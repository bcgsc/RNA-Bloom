/* 
 * Copyright (C) 2018-present BC Cancer Genome Sciences Centre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package rnabloom.util;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import rnabloom.io.ExtendedPafRecord;
import rnabloom.io.PafReader;
import rnabloom.io.PafRecord;
import rnabloom.olc.ComparableInterval;
import rnabloom.olc.Interval;

/**
 *
 * @author Ka Ming Nip
 */
public class PafUtils {
    private final static Pattern CIGAR_OP_PATTERN = Pattern.compile("(\\d+)([MIDNSHPX=])");
    
    public static boolean hasLargeOverlap(PafRecord r, int minOverlapMatches) {
        return (r.qEnd - r.qStart) >= minOverlapMatches &&
               (r.tEnd - r.tStart) >= minOverlapMatches;
    }
    
    public static boolean hasSimilarSizedOverlap(PafRecord r, int tolerance) {
        return Math.abs((r.qEnd - r.qStart) - (r.tEnd - r.tStart)) <= tolerance;
    }
    
    public static boolean hasGoodOverlap(PafRecord r, float minAlnId) {
        return r.numMatch / (float) r.blockLen >= minAlnId;
//        return r.numMatch / (float)(r.qEnd - r.qStart) >= minAlnId &&
//                r.numMatch / (float)(r.tEnd - r.tStart) >= minAlnId;
    }
            
    public static boolean hasAlignment(ExtendedPafRecord r) {
        return r.cigar != null && r.nm >= 0;
    }
    
    public static boolean hasGoodMatches(ExtendedPafRecord r, float minAlnId) {
        int totalCigarMatch = 0;
                
        Matcher m = CIGAR_OP_PATTERN.matcher(r.cigar);
        while (m.find()) {
            int opSize = Integer.parseInt(m.group(1));
            char op = m.group(2).charAt(0);
            if (op == 'M') {
                totalCigarMatch += opSize;
            }
        }
        
        float alnId = r.numMatch/(float)totalCigarMatch;
        
        return alnId >= minAlnId;
    }
    
    public static boolean hasGoodAlignment(ExtendedPafRecord r, int maxIndelSize, float minAlnId) {
        int totalCigarMatch = 0;
                
        Matcher m = CIGAR_OP_PATTERN.matcher(r.cigar);
        while (m.find()) {
            int opSize = Integer.parseInt(m.group(1));
            char op = m.group(2).charAt(0);
            switch (op) {
                case 'M':
                    totalCigarMatch += opSize;
                    break;
                case 'I':
                    if (opSize > maxIndelSize) {
                        return false;
                    }
                    break;
                case 'D':
                    if (opSize > maxIndelSize) {
                        return false;
                    }
                    break;
            }
        }
                
        return r.numMatch/(float)totalCigarMatch >= minAlnId && totalCigarMatch/(float)r.blockLen >= minAlnId;
    }
    
    public static boolean hasReverseComplementArtifact(PafRecord r, int maxEdgeClip) {
        if (r.qName.equals(r.tName) && r.reverseComplemented) {
            if (r.qStart <= maxEdgeClip || r.qLen - r.qEnd <= maxEdgeClip) {
                return true;
            }
        }
        
        return false;
    }
    
    public static enum CONTAIN_STATUS {QUERY, TARGET, BOTH, NEITHER};
    
    public static CONTAIN_STATUS getContained(PafRecord r, int maxEdgeClip) {
        
        boolean qContained = r.qStart <= maxEdgeClip && r.qLen - r.qEnd <= maxEdgeClip;
        boolean tContained = r.tStart <= maxEdgeClip && r.tLen - r.tEnd <= maxEdgeClip;
        
        if (qContained && tContained) {
            return CONTAIN_STATUS.BOTH;
        }
        else if (qContained) {
            return CONTAIN_STATUS.QUERY;
        }
        else if (tContained) {
            return CONTAIN_STATUS.TARGET;
        }
        
        return CONTAIN_STATUS.NEITHER;
    }
    
    public static boolean isQueryContained(PafRecord r, int maxEdgeClip) {
        return r.qStart <= maxEdgeClip && r.qLen - r.qEnd <= maxEdgeClip;
    }
    
    public static boolean isContainmentPafRecord(PafRecord r, int maxEdgeClip) {
        return ((r.qStart <= maxEdgeClip && r.qLen - r.qEnd <= maxEdgeClip) ||
            (r.tStart <= maxEdgeClip && r.tLen - r.tEnd <= maxEdgeClip));
    }
        
    public static boolean isForwardContainmentPafRecord(PafRecord r, int maxEdgeClip) {
        return !r.reverseComplemented && isContainmentPafRecord(r, maxEdgeClip);
    }
        
    public static boolean isDovetailPafRecord(PafRecord r, int maxEdgeClip) {
        if (r.reverseComplemented) {
            return (r.qEnd >= r.qLen - maxEdgeClip && r.tEnd >= r.tLen - maxEdgeClip && r.qStart > r.tLen - r.tEnd) ||
                    (r.tStart <= maxEdgeClip && r.qStart <= maxEdgeClip && r.qLen - r.qStart > r.tStart);
        }
        else {
            return (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip && r.qStart > r.tStart) ||
                    (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip && r.tStart > r.qStart);
        }
    }
    
    public static boolean isForwardDovetailPafRecord(PafRecord r, int maxEdgeClip) {
        if (!r.reverseComplemented) {
            if ((r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) ||
                    (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip)) {
                return true;
            }
        }
        
        return false;
    }
    
    private static String getBestPartner(HashMap<String, Interval> map, Interval source, int tolerance) {
        String bestPartner = null;
        int bestScore = 0;
        
        for (Entry<String, Interval> e : map.entrySet()) {
            Interval merged = ComparableInterval.merge(source, e.getValue());
            if (merged != null) {
                int score = (merged.end - source.end) + (source.start - merged.start);
                if (score > bestScore && score >= tolerance) {
                    bestScore = score;
                    bestPartner = e.getKey();
                }
            }
        }
        
        return bestPartner;
    }
    
    public static HashMap<String, Float> getLengthNormalizedReadCounts(String pafPath, Set<String> skipSet) throws IOException {
        HashMap<String, Float> counts = new HashMap<>();
        
        PafReader reader = new PafReader(pafPath);
        for (PafRecord r = new PafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if (!skipSet.contains(r.tName)) {
                Float c = counts.get(r.tName);
                if (c == null) {
                    c = r.numMatch/(float) r.tLen;
                }
                else {
                    c += r.numMatch/(float) r.tLen;
                }
                counts.put(r.tName, c);
            }
        }
        reader.close();
        
        return counts;
    }
    
    /**
     * 
     * @param pafPath input PAF file
     * @param targets set of target names
     * @param tolerance minimum extension required for bridging partners
     * @return map of targets and their read counts
     * @throws IOException 
     */
    public static HashMap<String, Float> getReadCounts(String pafPath, Set<String> targets, int tolerance) throws IOException {
        HashMap<String, Float> counts = new HashMap<>(targets.size());
        
        PafReader reader = new PafReader(pafPath);
        
        String prevName = null;
        int largestOverlapSize = 0;
        Interval largestOverlap = null;
        String bestTarget = null;
        HashMap<String, Interval> targetNameQueryIntervalMap = new HashMap<>();
        for (PafRecord r = new PafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if (targets.contains(r.tName)) {
                if (!r.qName.equals(prevName)) {
                    if (prevName != null) {
                        float increment = 1f;
                        targetNameQueryIntervalMap.remove(bestTarget);

                        String partner = getBestPartner(targetNameQueryIntervalMap, largestOverlap, tolerance);
                        if (partner != null) {
                            increment = 0.5f;
                            
                            // TODO: Store number of reads bridging two targets
                            if (counts.containsKey(partner)) {
                                float c = counts.get(partner);
                                counts.put(partner, c + increment);
                            }
                            else {
                                counts.put(partner, increment);
                            }
                        }
                        
                        if (counts.containsKey(bestTarget)) {
                            float c = counts.get(bestTarget);
                            counts.put(bestTarget, c + increment);
                        }
                        else {
                            counts.put(bestTarget, increment);
                        }
                        
                        targetNameQueryIntervalMap.clear();
                    }
                    
                    prevName = r.qName;
                    largestOverlapSize = 0;
                    bestTarget = null;
                }
                
                Interval i = new Interval(r.qStart, r.qEnd);
                targetNameQueryIntervalMap.put(r.tName, i);
                
                int overlapSize = r.qEnd - r.qStart;
                if (overlapSize > largestOverlapSize) {
                    largestOverlapSize = overlapSize;
                    bestTarget = r.tName;
                    largestOverlap = i;
                }
            }
        }
        
        reader.close();
        
        float increment = 1f;
        targetNameQueryIntervalMap.remove(bestTarget);

        String partner = getBestPartner(targetNameQueryIntervalMap, largestOverlap, tolerance);
        if (partner != null) {
            increment = 0.5f;
            
            if (counts.containsKey(partner)) {
                float c = counts.get(partner);
                counts.put(partner, c + increment);
            }
            else {
                counts.put(partner, increment);
            }
        }
        
        if (bestTarget != null) {
            if (counts.containsKey(bestTarget)) {
                float c = counts.get(bestTarget);
                counts.put(bestTarget, c + increment);
            }
            else {
                counts.put(bestTarget, increment);
            }
        }
        
        return counts;
    }
    
    public static void main(String[] args) {
        //debug
    }
}
