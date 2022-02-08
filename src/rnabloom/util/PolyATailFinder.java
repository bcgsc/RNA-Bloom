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
import static rnabloom.util.SeqUtils.reverseComplement;

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
    
    private final int seedLength;
    private final float minIdentity;
    private final int maxGap;
    private final int window;
    private final int pasSearchStart;
    private final int pasSearchEnd;
    
    public PolyATailFinder(int seedLength, float minIdentity, int maxGap, int window,
            int pasSearchStart, int pasSearchEnd) {
        this.seedLength = seedLength;
        this.minIdentity = minIdentity;
        this.maxGap = maxGap;
        this.window = window;
        this.pasSearchStart = pasSearchStart;
        this.pasSearchEnd = pasSearchEnd;
    }
    
    public static Pattern getPolyASignalPattern() {
        return Pattern.compile(String.join("|", POLY_A_SIGNALS), Pattern.CASE_INSENSITIVE);
    }
    
    public static Pattern getPolyASignalPatternRC() {
        return Pattern.compile(String.join("|", POLY_A_SIGNALS_REV_COMP), Pattern.CASE_INSENSITIVE);
    }
    
    public ArrayDeque<Integer> getPolyASignalPositions(String seq, int cleavagePos) {
        ArrayDeque<Integer> positions = new ArrayDeque<>();

        Matcher pasMatcher = pasPattern.matcher(seq);
        int start = cleavagePos + pasSearchStart;
        int end = cleavagePos + pasSearchEnd;
        pasMatcher.region(start, end);

        while (start < end && pasMatcher.find()) {
            start = pasMatcher.start();
            positions.add(start);
            pasMatcher.region(++start, end);
        }
        
        return positions;
    }

    public ArrayDeque<Integer> getPolyASignalPositionsRC(String seq, int cleavagePos) {
        ArrayDeque<Integer> positions = new ArrayDeque<>();

        Matcher pasMatcher = pasPatternRC.matcher(seq);
        int start = cleavagePos - pasSearchEnd;
        int end = cleavagePos - pasSearchStart;
        pasMatcher.region(start, end);

        while (start < end && pasMatcher.find()) {
            start = pasMatcher.start();
            positions.add(start);
            pasMatcher.region(++start, end);
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
        
        while (best != null && searchStart > 0) {
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
        String seq = "TCTGCTTGTCATCCACACAACACCGGGACTTAAGACAAATGGGACTGATGTCATCTTGAGCTCTTCATTTATTTTGACTGTGATTTATTTGGAGTGGAGGCATTGTTTCTAAGAAAAACATGTCATGTAGGTTGTCTAAAATAAAATGCGTTTTAAACTCAAAAAAAAAAAAAAAAA";
        int seedLength = 10;
        float minIdentity = 0.7f;
        int maxGap = 4;
        int window = 100;
        int pasSearchStart = -60;
        int pasSearchEnd = -5;
        
        PolyATailFinder finder = new PolyATailFinder(seedLength, minIdentity, maxGap, window, pasSearchStart, pasSearchEnd);
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
        
        String seqRC = reverseComplement(seq);
        itv = finder.findPolyTHead(seqRC);
        if (itv != null) {
            System.out.print(seqRC.substring(0, itv.start));
            System.out.print(seqRC.substring(itv.start, itv.end).toLowerCase());
            System.out.print(seqRC.substring(itv.end));
            System.out.println();
            System.out.flush();
            
            for (int pasPos : finder.getPolyASignalPositionsRC(seqRC, itv.end)) {
                System.out.println(Integer.toString(pasPos) + ' ' + seqRC.substring(pasPos, pasPos+6));
            }
        }
        
        System.out.flush();
    }
}
