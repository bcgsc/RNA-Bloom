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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import rnabloom.io.ExtendedPafRecord;
import rnabloom.io.PafReader;
import rnabloom.io.PafRecord;
import rnabloom.olc.Interval;
import static rnabloom.util.FileUtils.fileToStringCollection;
import static rnabloom.util.IntervalUtils.getOverlap;
import static rnabloom.util.IntervalUtils.merge;

/**
 *
 * @author Ka Ming Nip
 */
public class PafUtils {
    private final static Pattern CIGAR_OP_PATTERN = Pattern.compile("(\\d+)([MIDNSHPX=])");
    
    public static boolean hasLargeOverlap(PafRecord r, int minOverlapMatches) {
        return r.qEnd - r.qStart >= minOverlapMatches &&
               r.tEnd - r.tStart >= minOverlapMatches;
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

    public static boolean isTargetContained(PafRecord r, int maxEdgeClip) {
        return r.tStart <= maxEdgeClip && r.tLen - r.tEnd <= maxEdgeClip;
    }
    
    public static boolean isContainmentPafRecord(PafRecord r, int maxEdgeClip) {
        return isQueryContained(r, maxEdgeClip) || isTargetContained(r, maxEdgeClip);
    }
        
    public static boolean isForwardContainmentPafRecord(PafRecord r, int maxEdgeClip) {
        return !r.reverseComplemented && isContainmentPafRecord(r, maxEdgeClip);
    }
    
    public static boolean isQueryEdgeSink(PafRecord r, int maxEdgeClip) {
        int qHead = r.qStart;
        int qTail = r.qLen - r.qEnd;
        int tHead = r.tStart;
        int tTail = r.tLen - r.tEnd;
        
        if (r.reverseComplemented) {
            if (qTail > tHead && tTail > qHead && qHead <= maxEdgeClip && tHead <= maxEdgeClip) {
                return true;
            }
        }
        else {
            if (tHead > qHead && tTail < qTail && qHead <= maxEdgeClip && tTail <= maxEdgeClip) {
                return true;
            }
        }
        
        return false;
    }
    
    public static boolean isDovetailPafRecord(PafRecord r, int maxEdgeClip) {
        int qHead = r.qStart;
        int qTail = r.qLen - r.qEnd;
        int tHead = r.tStart;
        int tTail = r.tLen - r.tEnd;

        if (r.reverseComplemented) {
            return ((qHead > tTail && tHead > qTail && tTail <= maxEdgeClip && qTail <= maxEdgeClip) ||
                    (qTail > tHead && tTail > qHead && qHead <= maxEdgeClip && tHead <= maxEdgeClip));
        }
        else {
            return ((qHead > tHead && qTail < tTail && tHead <= maxEdgeClip && qTail <= maxEdgeClip) ||
                    (tHead > qHead && tTail < qTail && qHead <= maxEdgeClip && tTail <= maxEdgeClip));
        }
    }
    
    public static boolean isForwardDovetailPafRecord(PafRecord r, int maxEdgeClip) {
        if (!r.reverseComplemented) {
            int qHead = r.qStart;
            int qTail = r.qLen - r.qEnd;
            int tHead = r.tStart;
            int tTail = r.tLen - r.tEnd;
            
            return (qHead > tHead && qTail < tTail && tHead <= maxEdgeClip && qTail <= maxEdgeClip) ||
                    (tHead > qHead && tTail < qTail && qHead <= maxEdgeClip && tTail <= maxEdgeClip);
        }
        
        return false;
    }
    
    private static String getBestPartner(HashMap<String, Interval> map, Interval source, int tolerance) {
        String bestPartner = null;
        int bestScore = 0;
        
        for (Entry<String, Interval> e : map.entrySet()) {
            Interval merged = merge(source, e.getValue());
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
    
    private static class Overlap implements Comparable<Overlap> {
        String tName;
        int numMatch, tLen;
        int qStart, qEnd, tStart, tEnd;
        
        public Overlap(PafRecord r) {
            this.tName = r.tName;
            this.numMatch = r.numMatch;
            this.tLen = r.tLen;
            this.qStart = r.qStart;
            this.qEnd = r.qEnd;
            this.tStart = r.tStart;
            this.tEnd = r.tEnd;
        }
        
        @Override
        public int compareTo(Overlap other) {
            return (other.qEnd - other.qStart) - (this.qEnd - this.qStart);
        }
    }
    
    private static boolean isOverlapContained(Overlap q, Collection<Overlap> list) {
        for (Overlap m : list) {
            if (q.qStart >= m.qStart && q.qEnd <= m.qEnd && q.numMatch < m.numMatch * 0.95f) {
                return true;
            }
        }
        return false;
    }
    
    private static Overlap getOverlapContainer(Overlap q, Collection<Overlap> list) {
        for (Overlap m : list) {
            if (q.qStart >= m.qStart && q.qEnd <= m.qEnd) {
                return m;
            }
        }
        return null;
    }
    
    private static Overlap getOverlapContainer(Overlap q, Collection<Overlap> list, float maxProportion) {
        int maxOverlapLen = 0;
        Overlap container = null;
        for (Overlap other : list) {
            int overlapLen = getOverlap(q.qStart, q.qEnd, other.qStart, other.qEnd);
            if (overlapLen > maxOverlapLen) {
                maxOverlapLen = overlapLen;
                container = other;
            }
        }
        
        if (maxOverlapLen >= maxProportion * (q.qEnd - q.qStart)) {
            return container;
        }
        
        return null;
    }
    
    private static void updateNormalizedCounts(ArrayList<Overlap> targets, HashMap<String, Float> counts) {
        if (!targets.isEmpty()) {
            if (targets.size() == 1) {
                Overlap t = targets.get(0);
                Float c = counts.get(t.tName);
                if (c == null) {
                    c = (t.tEnd - t.tStart)/(float) t.tLen;
                }
                else {
                    c += (t.tEnd - t.tStart)/(float) t.tLen;
                }
                counts.put(t.tName, c);
            }
            else {
                Collections.sort(targets);

                ArrayList<Overlap> targetsKept = new ArrayList<>();
                HashMap<Overlap, ArrayList<Overlap>> multiTargets = new HashMap<>();
                for (Overlap m : targets) {
                    Overlap c = getOverlapContainer(m, targetsKept, 0.95f);
                    if (c == null) {
                        // region not contained
                        targetsKept.add(m);
                    }
                    else if (m.qEnd - m.qStart >= (c.qEnd - c.qStart) * 0.95f) {
                        // region multimaps
                        ArrayList<Overlap> multimaps = multiTargets.get(c);
                        if (multimaps == null) {
                            multimaps = new ArrayList<>();
                            multiTargets.put(c, multimaps);
                        }
                        multimaps.add(m);
                    }
                }

                for (Overlap t : targetsKept) {
                    ArrayList<Overlap> multimaps = multiTargets.get(t);
                    if (multimaps != null) {
                        // fractional assignment of multimapped region
                        float fraction = 1f/(multimaps.size() + 1);

                        Float c = counts.get(t.tName);
                        if (c == null) {
                            c = (t.tEnd - t.tStart)/(float) t.tLen * fraction;
                        }
                        else {
                            c += (t.tEnd - t.tStart)/(float) t.tLen * fraction;
                        }
                        counts.put(t.tName, c);

                        for (Overlap mm : multimaps) {
                            c = counts.get(mm.tName);
                            if (c == null) {
                                c = (mm.tEnd - mm.tStart)/(float) mm.tLen * fraction;
                            }
                            else {
                                c += (mm.tEnd - mm.tStart)/(float) mm.tLen * fraction;
                            }
                            counts.put(mm.tName, c);
                        }
                    }
                    else {
                        Float c = counts.get(t.tName);
                        if (c == null) {
                            c = (t.tEnd - t.tStart)/(float) t.tLen;
                        }
                        else {
                            c += (t.tEnd - t.tStart)/(float) t.tLen;
                        }
                        counts.put(t.tName, c);
                    }
                }
            }
        }
    }
    
    public static HashMap<String, Float> getLengthNormalizedReadCounts(String pafPath, Set<String> skipSet) throws IOException {
        HashMap<String, Float> counts = new HashMap<>();
                
        PafReader reader = new PafReader(pafPath);
        String prevName = null;
        ArrayList<Overlap> targets = new ArrayList<>();
        for (PafRecord r = new PafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if (!r.qName.equals(prevName)) {
                updateNormalizedCounts(targets, counts);
                targets = new ArrayList<>();
                prevName = r.qName;
            }
            
            if (!skipSet.contains(r.tName)) {
                targets.add(new Overlap(r));
            }
        }
        reader.close();
        
        // process the last batch of records
        updateNormalizedCounts(targets, counts);
        
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
