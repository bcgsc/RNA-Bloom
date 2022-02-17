/* 
 * Copyright (C) 2022-present BC Cancer Genome Sciences Centre
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

import java.util.ArrayDeque;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import rnabloom.olc.Interval;

/**
 *
 * @author Ka Ming Nip
 */
public class PolyATailFinder {    
    public static final String[] POLY_A_SIGNALS = {"AATAAA", "ATTAAA", "AGTAAA", "TATAAA",
                                                   "CATAAA", "GATAAA", "AATATA", "AATACA",
                                                   "AATAGA", "AAAAAG", "ACTAAA", "AAGAAA",
                                                   "AATGAA", "TTTAAA", "AAAACA", "GGGGCT",
                                                   "AATAAT", "AACAAA", "ATTACA", "ATTATA",
                                                   "AACAAG", "AATAAG", "TTTTTT"}; // PMID: 27382025

    public static final String[] POLY_A_SIGNALS_REV_COMP = {"TTTATT", "TTTAAT", "TTTACT", "TTTATA",
                                                   "TTTATG", "TTTATC", "TATATT", "TGTATT",
                                                   "TCTATT", "CTTTTT", "TTTAGT", "TTTCTT",
                                                   "TTCATT", "TTTAAA", "TGTTTT", "AGCCCC",
                                                   "ATTATT", "TTTGTT", "TGTAAT", "TATAAT",
                                                   "CTTGTT", "CTTATT", "AAAAAA"};
    
    private final Pattern pasPattern = getPolyASignalPattern();
    private final Pattern pasPatternRC = getPolyASignalPatternRC();
    
    private int seedLength = 10;
    private float minIdentity = 0.7f;
    private int maxGap = 4;
    private int window = 100; // search window for polyA seed
    private int pasSearchStart = 60; // distance upstream of cleavage site position
    private int pasSearchEnd = 5;
    
    public static enum Profile {ILLUMINA, ONT};
    
    public PolyATailFinder() {
        
    }
            
    public PolyATailFinder(int seedLength, float minIdentity, int maxGap, int window,
            int pasSearchStart, int pasSearchEnd) {
        this.seedLength = seedLength;
        this.minIdentity = minIdentity;
        this.maxGap = maxGap;
        this.window = window;
        this.pasSearchStart = pasSearchStart;
        this.pasSearchEnd = pasSearchEnd;
    }
    
    public void setProfile(Profile p) {
        switch(p) {
            case ILLUMINA:
                seedLength = 4;
                minIdentity = 0.9f;
                maxGap = 1;
                window = seedLength+maxGap;
                pasSearchStart = 60;
                pasSearchEnd = 5;
                break;
            case ONT:
                seedLength = 12;
                minIdentity = 0.9f;
                maxGap = 4;
                window = 30;
                pasSearchStart = 60;
                pasSearchEnd = 5;
                break;
        }
        window = Math.max(seedLength, window);
    }
    
    public void setSeedLength(int seedLength) {
        this.seedLength = seedLength;
    }
    
    public int getSeedLength() {
        return seedLength;
    }
    
    public void setWindow(int window) {
        this.window = Math.max(seedLength, window);
    }
    
    public void setMinIdentity(float minIdentity) {
        this.minIdentity = minIdentity;
    }
    
    public void setMaxGap(int maxGap) {
        this.maxGap = maxGap;
    }
    
    public void setPasSearchRange(int start, int end) {
        this.pasSearchStart = start;
        this.pasSearchEnd = end;
    }
    
    public static Pattern getPolyASignalPattern() {
        return Pattern.compile(String.join("|", POLY_A_SIGNALS), Pattern.CASE_INSENSITIVE);
    }
    
    public static Pattern getPolyASignalPatternRC() {
        return Pattern.compile(String.join("|", POLY_A_SIGNALS_REV_COMP), Pattern.CASE_INSENSITIVE);
    }
    
    public boolean hasPolyASignal(String seq, int cleavagePos) {
        Matcher pasMatcher = pasPattern.matcher(seq);
        int start = Math.max(0, cleavagePos - pasSearchStart);
        int end = Math.max(0, cleavagePos - pasSearchEnd);
        if (start < end) {
            pasMatcher.region(start, end);
            return pasMatcher.find();
        }
        return false;
    }
    
    public boolean hasPolyASignalRC(String seq, int cleavagePos) {
        int seqLength = seq.length();
        Matcher pasMatcher = pasPatternRC.matcher(seq);
        int start = Math.min(cleavagePos + pasSearchEnd, seqLength);
        int end = Math.min(cleavagePos + pasSearchStart, seqLength);
        if (start < end) {
            pasMatcher.region(start, end);
            return pasMatcher.find();
        }
        return false;
    }
    
    public ArrayDeque<Integer> getPolyASignalPositions(String seq, int cleavagePos) {
        ArrayDeque<Integer> positions = new ArrayDeque<>();

        Matcher pasMatcher = pasPattern.matcher(seq);
        int start = Math.max(0, cleavagePos - pasSearchStart);
        int end = Math.max(0, cleavagePos - pasSearchEnd);
        if (start < end) {
            pasMatcher.region(start, end);

            while (start < end && pasMatcher.find()) {
                start = pasMatcher.start();
                positions.add(start);
                pasMatcher.region(++start, end);
            }
        }
        
        return positions;
    }

    public ArrayDeque<Integer> getPolyASignalPositionsRC(String seq, int cleavagePos) {
        ArrayDeque<Integer> positions = new ArrayDeque<>();

        int seqLength = seq.length();
        Matcher pasMatcher = pasPatternRC.matcher(seq);
        int start = Math.min(cleavagePos + pasSearchEnd, seqLength);
        int end = Math.min(cleavagePos + pasSearchStart, seqLength);
        if (start < end) {
            pasMatcher.region(start, end);

            while (start < end && pasMatcher.find()) {
                start = pasMatcher.start();
                positions.add(start);
                pasMatcher.region(++start, end);
            }
        }
        
        return positions;
    }
    
    private static int countA(String seq, int start, int end) {
        int count = 0;
        for (int i=start; i<end; ++i) {
            char c = seq.charAt(i);
            if (c == 'A' || c == 'a') {
                ++count;
            }
        }
        return count;
    }
    
    private static float percentA(String seq, int start, int end) {
        return (float)(countA(seq, start, end))/(end - start);
    }
    
    public Interval findPolyASeed(String seq, int searchStart, int searchEnd) {
        if (searchStart < searchEnd && searchStart >=0 && searchEnd >= 0 && searchEnd-searchStart >= seedLength) {
            int numA = 0;
            for (int i=searchEnd-1; i>=searchEnd-seedLength; --i) {
                char c = seq.charAt(i);
                if (c == 'A' || c == 'a') {
                    ++numA;
                }
            }
            
            float id = numA/(float)seedLength;
            Interval bestRegion = null;
            
            if (id >= minIdentity) {
                bestRegion = new Interval(searchEnd-seedLength, searchEnd);
            }
            
            for (int i=searchEnd-seedLength-1; i>=searchStart; --i) {
                if (numA > 0) {
                    char out = seq.charAt(i+seedLength);
                    if (out == 'A' || out == 'a') {
                        --numA;
                    }
                }
                
                char c = seq.charAt(i);
                if (c == 'A' || c == 'a') {
                    id = ++numA/(float)seedLength;
                    if (bestRegion == null) {
                        if (id >= minIdentity) {
                            bestRegion = new Interval(i, i+seedLength);
                        }
                    }
                    else {
                        if (id >= minIdentity) {
                            bestRegion.start = i;
                        }
                        else {
                            break;
                        }
                    }
                }
                else if (bestRegion != null && numA/(float)seedLength < minIdentity){
                    break;
                }
            }
            
            if (bestRegion != null) {
                // refine region
                while (bestRegion.end - bestRegion.start > seedLength) {
                    char c = seq.charAt(bestRegion.end-1);
                    if (c == 'A' || c == 'a') {
                        break;
                    }
                    else {
                        --bestRegion.end;
                    }
                }
                
                if (bestRegion.start == searchStart) {
                    while (bestRegion.start > 0) {
                        char c = seq.charAt(bestRegion.start);
                        if (c == 'A' || c == 'a') {
                            --bestRegion.start;
                        }
                        else {
                            break;
                        }
                    }
                }
            }
            
            return bestRegion;
        }
        
        return null;
    }
    
    public Interval findPolyAPerfectSeed(String seq, int searchStart, int searchEnd) {
        if (searchStart < searchEnd && searchStart >=0 && searchEnd >= 0) {
            int memStart = -1;
            int numA = 0;
            for (int i=searchEnd-1; i>=searchStart; --i) {
                char c = seq.charAt(i);
                if (c == 'A' || c == 'a') {
                    if (++numA >= seedLength) {
                        memStart = i;
                    }
                }
                else if (numA < seedLength)  {
                    numA = 0;
                }
                else {
                    break;
                }
            }

            if (memStart > 0) {
                if (memStart == searchStart) {
                    for (int i=searchStart-1; i>=0; --i) {
                        char c = seq.charAt(i);
                        if (c == 'A' || c == 'a') {
                            ++numA;
                            memStart = i;
                        }
                        else {
                            break;
                        }
                    }
                }

                return new Interval(memStart, memStart + numA);
            }
        }
        
        return null;
    }
    
    public Interval findPolyATail(String seq) {
        int seqLen = seq.length();
        int searchEnd = seqLen;
        int searchStart = Math.max(0, searchEnd - window);
        Interval best = findPolyASeed(seq, searchStart, searchEnd);
        
        while (best != null && searchStart > 0) {
            searchEnd = best.start;
            searchStart = Math.max(0, searchEnd - window);
            Interval prev = findPolyASeed(seq, searchStart, searchEnd);
            if (prev != null && (prev.end + maxGap >= best.start || percentA(seq, prev.end, best.start) >= minIdentity)) {
                best.start = prev.start;
            }
            else {
                break;
            }
        }
        
        return best;
    }
    
    private static int countT(String seq, int start, int end) {
        int count = 0;
        for (int i=start; i<end; ++i) {
            char c = seq.charAt(i);
            if (c == 'T' || c == 't') {
                ++count;
            }
        }
        return count;
    }
    
    private static float percentT(String seq, int start, int end) {
        return (float)(countT(seq, start, end))/(end - start);
    }
    
    public Interval findPolyTSeed(String seq, int searchStart, int searchEnd) {
        if (searchStart < searchEnd && searchStart >=0 && searchEnd >= 0 && searchEnd-searchStart >= seedLength) {
            int numT = 0;
            for (int i=searchStart; i<searchStart+seedLength; ++i) {
                char c = seq.charAt(i);
                if (c == 'T' || c == 't') {
                    ++numT;
                }
            }
            
            float id = numT/(float)seedLength;
            Interval bestRegion = null;
            
            if (id >= minIdentity) {
                bestRegion = new Interval(searchStart, searchStart+seedLength);
            }
            
            for (int i=searchStart+seedLength; i<searchEnd; ++i) {
                if (numT > 0) {
                    char out = seq.charAt(i-seedLength);
                    if (out == 'T' || out == 't') {
                        --numT;
                    }
                }
                
                char c = seq.charAt(i);
                if (c == 'T' || c == 't') {
                    id = ++numT/(float)seedLength;
                    if (bestRegion == null) {
                        if (id >= minIdentity) {
                            bestRegion = new Interval(i-seedLength, i+1);
                        }
                    }
                    else {
                        if (id >= minIdentity) {
                            bestRegion.end = i+1;
                        }
                        else {
                            break;
                        }
                    }
                }
                else if (bestRegion != null && numT/(float)seedLength < minIdentity){
                    break;
                }
            }
            
            if (bestRegion != null) {
                // refine region
                while (bestRegion.end - bestRegion.start > seedLength) {
                    char c = seq.charAt(bestRegion.start);
                    if (c == 'T' || c == 't') {
                        break;
                    }
                    else {
                        ++bestRegion.start;
                    }
                }
                
                if (bestRegion.end == searchEnd) {
                    int seqLen = seq.length();
                    while (bestRegion.end < seqLen) {
                        char c = seq.charAt(bestRegion.end);
                        if (c == 'T' || c == 't') {
                            ++bestRegion.end;
                        }
                        else {
                            break;
                        }
                    }
                }
            }
            
            return bestRegion;
        }
        
        return null;
    }
    
    public Interval findPolyTPerfectSeed(String seq, int searchStart, int searchEnd) {
        if (searchStart < searchEnd && searchStart >=0 && searchEnd >= 0) {
            int seqLen = seq.length();
            
            int memEnd = -1;
            int numT = 0;
            for (int i=searchStart; i<searchEnd; ++i) {
                char c = seq.charAt(i);
                if (c == 'T' || c == 't') {
                    if (++numT >= seedLength) {
                        memEnd = i;
                    }
                }
                else if (numT < seedLength)  {
                    numT = 0;
                }
                else {
                    break;
                }
            }

            if (memEnd > 0) {
                if (memEnd == searchEnd-1) {
                    for (int i=searchEnd; i<seqLen; ++i) {
                        char c = seq.charAt(i);
                        if (c == 'T' || c == 't') {
                            ++numT;
                            memEnd = i;
                        }
                        else {
                            break;
                        }
                    }
                }

                return new Interval(memEnd-numT+1, memEnd+1);
            }
        }
        
        return null;
    }
    
    public Interval findPolyTHead(String seq) {
        int seqLen = seq.length();
        int searchStart = 0;
        int searchEnd = Math.min(seqLen, window);
        Interval best = findPolyTSeed(seq, searchStart, searchEnd);
        
        while (best != null && searchEnd < seqLen) {
            searchStart = best.end;
            searchEnd = Math.min(seqLen, best.end + window);
            Interval next = findPolyTSeed(seq, searchStart, searchEnd);
            if (next != null && (best.end + maxGap >= next.start || percentT(seq, best.end, next.start) >= minIdentity)) {
                best.end = next.end;
            }
            else {
                break;
            }
        }
        
        return best;
    }
    
    public static void main(String[] args) {
        //debug
        String seq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCTGGAATTGAATGTCTTTATTAAATAAACGAGTAAATGGTAGCACAAATCACCATCAATATTTTTGGAAGGATTGGGGACAAGATGTCGAGTCAGAATATAATTGTTCATTTCAGGGTCTCAATGTAGCTGAAGAACTGTGCCCACTGATCGAGTATTAGGTATTGCAAATGCAGGAGGTAAGGCTAAGAAATAGGACTTACGCCGTTCAGAAGAGATTGAATTGAAACCTTAAAAACTATCATAATAGTAGGAATGCATGTTAAGATTTGATAACTTTCTTTAACTAGAGTTTTCAACCCACAGTTAGGAGCAAAGTTGTAAAGTGAGTAGGTGTGAAGAAGGACACTCTTTTGAAAGAAATTAACTACTTCTAAAATGATTTCAATTATTTTCCCTATTTTCATTTTCCATAATATTTCCTCCTGTTTTATACCTCTAACTACCCTCTGATTTCTCTGAGGAAAAAAAAAAGAATATAAGAGCAGGATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        System.out.println(seq);
        PolyATailFinder finder = new PolyATailFinder();
        finder.setProfile(Profile.ONT);
        Interval itv = finder.findPolyATail(seq);
        if (itv != null) {
            System.out.print(seq.substring(0, itv.start));
            System.out.print(seq.substring(itv.start, itv.end).toLowerCase());
            System.out.print(seq.substring(itv.end));
            System.out.println();
            System.out.flush();
            
            for (int pasPos : finder.getPolyASignalPositions(seq, itv.start)) {
                System.out.println(Integer.toString(pasPos) + ' ' + seq.substring(pasPos, pasPos+6));
            }
        }
        
        //String seqRC = reverseComplement(seq);
        itv = finder.findPolyTHead(seq);
        if (itv != null) {
            System.out.print(seq.substring(0, itv.start));
            System.out.print(seq.substring(itv.start, itv.end).toLowerCase());
            System.out.print(seq.substring(itv.end));
            System.out.println();
            System.out.flush();
            
            for (int pasPos : finder.getPolyASignalPositionsRC(seq, itv.end)) {
                System.out.println(Integer.toString(pasPos) + ' ' + seq.substring(pasPos, pasPos+6));
            }
        }
        
        System.out.flush();
    }
}
