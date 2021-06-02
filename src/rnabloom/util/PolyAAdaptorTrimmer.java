/* 
 * Copyright (C) 2021-present BC Cancer Genome Sciences Centre
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

/**
 *
 * @author Ka Ming Nip
 */
public class PolyAAdaptorTrimmer {
    private final Pattern polyATailPattern;
    private final Pattern polyTHeadPattern;
    private final int minPolyALen;
    
    public PolyAAdaptorTrimmer(int minPolyALen) {
        this(minPolyALen, false);
    }
    
    public PolyAAdaptorTrimmer(int minPolyALen, boolean caseInsensitive) {
        this.minPolyALen = minPolyALen;
        
        if (caseInsensitive) {
            polyATailPattern = Pattern.compile(".*A{" + minPolyALen + ",}", Pattern.CASE_INSENSITIVE);
            polyTHeadPattern = Pattern.compile("T{" + minPolyALen + ",}", Pattern.CASE_INSENSITIVE);
        }
        else {
            polyATailPattern = Pattern.compile(".*A{" + minPolyALen + ",}");
            polyTHeadPattern = Pattern.compile("T{" + minPolyALen + ",}");
        }
    }
    
    public String chompHeadAdaptor(String seq, int searchRange) {
        int seqLen = seq.length();
        
        if (seqLen >= minPolyALen) {
            Matcher m = polyTHeadPattern.matcher(seq);
            if (searchRange < seqLen) {
                m = m.region(0, searchRange);
            }

            if (m.find()) {
                int start = m.start();
                if (start > 0) {
                    return seq.substring(start);
                }
            }
        }
        
        return seq;
    }
    
    public String chompTailAdaptor(String seq, int searchRange) {
        int seqLen = seq.length();
        
        if (seqLen >= minPolyALen) {
            Matcher m = polyATailPattern.matcher(seq);
            int start = seqLen - searchRange;
            if (start > 0) {
                m = m.region(start, seqLen);
            }

            if (m.find()) {
                int end = m.end();
                if (end < seqLen) {
                    return seq.substring(0, end);
                }
            }
        }
        
        return seq;
    }
    
    public static void main(String[] args) {
        //debug
    }
}
