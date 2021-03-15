/* 
 * Copyright (C) 2021 BC Cancer Genome Sciences Centre
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

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import rnabloom.io.ExtendedPafRecord;
import rnabloom.io.PafRecord;

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
    
    public static boolean hasGoodOverlap(PafRecord r, float minAlnId) {
        return r.numMatch / (float) r.blockLen >= minAlnId;
//        return r.numMatch / (float)(r.qEnd - r.qStart) >= minAlnId &&
//                r.numMatch / (float)(r.tEnd - r.tStart) >= minAlnId;
    }
        
    public static boolean hasAlignment(ExtendedPafRecord r) {
        return r.cigar != null && r.nm >= 0;
    }
    
    public static boolean hasGoodAlignment(ExtendedPafRecord r, int maxIndelSize, float minAlnId) {
        int numMatch = 0;
        //int numDel = 0;
        //int numIns = 0;
                
        Matcher m = CIGAR_OP_PATTERN.matcher(r.cigar);
        while (m.find()) {
            int opSize = Integer.parseInt(m.group(1));
            char op = m.group(2).charAt(0);
            switch (op) {
                case 'M':
                    numMatch += opSize;
                    break;
                case 'I':
                    if (opSize > maxIndelSize) {
                        return false;
                    }
                    //numIns += opSize;
                    break;
                case 'D':
                    if (opSize > maxIndelSize) {
                        return false;
                    }
                    //numDel += opSize;
                    break;
            }
        }
        
        //float alnId = (numMatch - record.nm)/(float)(numMatch + numDel + numIns);
        float alnId = (numMatch - r.nm)/(float)r.blockLen;
        
        return alnId >= minAlnId;
    }
    
    public static boolean hasReverseComplementArtifact(PafRecord r, int maxEdgeClip) {
        if (r.qName.equals(r.tName) && r.reverseComplemented) {
            if (r.qStart <= maxEdgeClip || r.qLen - r.qEnd <= maxEdgeClip) {
                return true;
            }
        }
        
        return false;
    }
    
    public static boolean isContainmentPafRecord(PafRecord r, int maxEdgeClip) {
        return ((r.qStart <= maxEdgeClip && r.qLen - r.qEnd <= maxEdgeClip) ||
            (r.tStart <= maxEdgeClip && r.tLen - r.tEnd <= maxEdgeClip));
    }
        
    public static boolean isStrandedContainmentPafRecord(PafRecord r, int maxEdgeClip) {
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
    
    public static boolean isStrandedDovetailPafRecord(PafRecord r, int maxEdgeClip) {
        if (!r.reverseComplemented) {
            if ((r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) ||
                    (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip)) {
                return true;
            }
        }
        
        return false;
    }
}
