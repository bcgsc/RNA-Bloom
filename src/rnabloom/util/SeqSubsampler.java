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

import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.hash.CanonicalHashFunction;
import rnabloom.bloom.hash.HashFunction;
import rnabloom.bloom.hash.MinimizerHashIterator;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import static rnabloom.util.Common.convertToRoundedPercent;
import static rnabloom.util.SeqUtils.chompPolyATail;
import static rnabloom.util.SeqUtils.chompPolyTHead;
import static rnabloom.util.SeqUtils.compressHomoPolymers;

/**
 *
 * @author Ka Ming Nip
 */
public class SeqSubsampler {
    
    public static void minimizerBased(String inFasta, String outFasta,
            long bfSize, int k, int w, int numHash, boolean stranded, 
            int maxNonMatchingChainLength, float minMatchingProportion,
            int maxMultiplicity) throws IOException {        
        // read all sequences
        ArrayList<BitSequence> seqs = new ArrayList<>();
        FastaReader fr = new FastaReader(inFasta);
        while(fr.hasNext()) {
            seqs.add(new BitSequence(fr.next()));
        }
        fr.close();
        int numSeq = seqs.size();
        
        // sort from longest to shortest
        Collections.sort(seqs);
        
        HashFunction h = stranded ? new HashFunction(k) : new CanonicalHashFunction(k);
        
        MinimizerHashIterator itr = new MinimizerHashIterator(k, w, h.getHashIterator(1));
        CountingBloomFilter bf = new CountingBloomFilter(bfSize, numHash, h);
        
        FastaWriter fw  = new FastaWriter(outFasta, false);
        int seqID = 0;
        long[] hVals = new long[numHash];
        for (BitSequence s : seqs) {
            String seq = s.toString();
            if (itr.start(seq)) {
                int numMinimizers = 0;
                int numMinimizersSeen = 0;
                int consecutiveMissing = 0;
                int maxConsecutiveMissing = 0;
                
                if (itr.hasNext()) {
                    long prev = itr.next();
                    ++numMinimizers;
                    itr.getMultipleHashValues(prev, hVals);
                    if (bf.incrementAndGet(hVals) > maxMultiplicity) {
                        ++numMinimizersSeen;
                    }
                    
                    while (itr.hasNext()) {
                        long mm = itr.next();
                        if (mm != prev) {
                            ++numMinimizers;
                            itr.getMultipleHashValues(mm, hVals);
                            if (bf.incrementAndGet(hVals) > maxMultiplicity) {
                                ++numMinimizersSeen;
                                consecutiveMissing = 0;
                            }
                            else {
                                maxConsecutiveMissing = Math.max(maxConsecutiveMissing, ++consecutiveMissing);
                            }
                        }
                        prev = mm;
                    }
                }
                
                if (maxConsecutiveMissing > maxNonMatchingChainLength || (float)numMinimizersSeen/(float)numMinimizers < minMatchingProportion) {
                    fw.write(Integer.toString(++seqID), seq);
                }
            }
            else {
                fw.write(Integer.toString(++seqID), seq);
            }
        }
        fw.close();
        
        System.out.println("Bloom filter FPR:\t" + convertToRoundedPercent(bf.getFPR()) + " %");
        
        bf.destroy();
        
        System.out.println("before: " + NumberFormat.getInstance().format(numSeq) + 
                            "\tafter: " + NumberFormat.getInstance().format(seqID) +
                            " (" + convertToRoundedPercent(seqID/(float)numSeq) + " %)");

    }
    
    public static void kmerBased(ArrayList<? extends BitSequence> seqs, String outSubsampleFasta,
            long bfSize, int k, int numHash, boolean stranded,
            int maxMultiplicity, int maxEdgeClip, boolean verbose) throws IOException {
        int numSeq = seqs.size();
        
        HashFunction h = stranded ? new HashFunction(k) : new CanonicalHashFunction(k);
        
        CountingBloomFilter cbf = new CountingBloomFilter(bfSize, numHash, h);
        NTHashIterator itr = h.getHashIterator(numHash, k);
        long[] hVals = itr.hVals;
        
        FastaWriter writer  = new FastaWriter(outSubsampleFasta, false);
                
        int numSubsample = 0;

        for (BitSequence s : seqs) {                    
            String seq = s.toString();
            
            int seqLen = seq.length();
            boolean tooShort = seqLen < 2 * maxEdgeClip + k;
            int maxPos = seqLen - maxEdgeClip;
            int maxMissingChainLen = 0;
            int missingChainLen = 0;
            
            if (itr.start(seq)) {
                if (tooShort) {
                    while (itr.hasNext()) {
                        itr.next();
                        if (cbf.incrementAndGet(hVals) <= maxMultiplicity) {
                            ++missingChainLen;
                        }
                        else {
                            if (missingChainLen > maxMissingChainLen) {
                                maxMissingChainLen = missingChainLen;
                            }
                            missingChainLen = 0;
                        }
                    }
                }
                else {
                    // head region
                    for (int pos=0; pos<maxEdgeClip; ++pos) {
                        itr.next();
                        cbf.increment(hVals);
                    }
                    
                    int effRangeEndIndex = maxPos-maxEdgeClip;
                    for (int pos=maxEdgeClip; pos<effRangeEndIndex; ++pos) {
                        itr.next();
                        if (cbf.incrementAndGet(hVals) <= maxMultiplicity) {
                            ++missingChainLen;
                        }
                        else {
                            if (missingChainLen > maxMissingChainLen) {
                                maxMissingChainLen = missingChainLen;
                            }
                            missingChainLen = 0;
                        }
                    }
                    
                    // tail region
                    for (int pos=effRangeEndIndex; pos<maxPos; ++pos) {
                        itr.next();
                        cbf.increment(hVals);
                    }
                }
                
                if (missingChainLen > maxMissingChainLen) {
                    maxMissingChainLen = missingChainLen;
                }
            }
            
            if (maxMissingChainLen >= k) {
                ++numSubsample;
                writer.write("s" + numSubsample, seq);
            }
        }
        
        writer.close();
        
        float fpr = cbf.getFPR();
        cbf.destroy();
        
        if (verbose) {
            System.out.println("Bloom filter FPR:\t" + convertToRoundedPercent(fpr) + " %");
            System.out.println("before: " + NumberFormat.getInstance().format(numSeq) + 
                                "\tafter: " + NumberFormat.getInstance().format(numSubsample) +
                                " (" + convertToRoundedPercent(numSubsample/(float)numSeq) + " %)");
        }
    }
    
    public static void minmalSet(ArrayList<? extends BitSequence> seqs, String outFasta,
            long bfSize, int k, int numHash, boolean stranded, boolean useHpcKmers,
            int windowSize, int minMatchingWindows, int minSeqLen) throws IOException, InterruptedException {
        int numSeq = seqs.size();
        
        FastaWriter writer = new FastaWriter(outFasta, false);
        
        HashFunction h = stranded ? new HashFunction(k) : new CanonicalHashFunction(k);
        
        BloomFilter bf = new BloomFilter(bfSize, numHash, h);
        NTHashIterator itr = h.getHashIterator(numHash, k);
        long[] hVals = itr.hVals;
        int id = 0;
        
        for (BitSequence s : seqs) {
            if (s != null) {
                String seq = s.toString();
                
                // chomp poly A tail and poly T head
                seq = chompPolyATail(seq, 50, 0.85f);
                if (!stranded) {
                    seq = chompPolyTHead(seq, 50, 0.85f);
                }
                
                if (seq.length() >= minSeqLen) {
                    String hpc = useHpcKmers ? compressHomoPolymers(seq) : seq;

                    if (itr.start(hpc)) {
    //                    int numKmers = hpc.length() - k + 1;
    //                    int numWindows = numKmers/windowSize;
    //                    if (numKmers % windowSize > 0) {
    //                        ++numWindows;
    //                    }

                        int numWindowsSeen = 0;
                        int windowIndex = 0;
                        boolean windowStatus = false;
                        while (itr.hasNext()) {
                            itr.next();

                            if (itr.getPos()/windowSize > windowIndex) {
                                ++windowIndex;
                                windowStatus = false;
                            }

                            if (bf.lookup(hVals)) {
                                if (!windowStatus) {
                                    ++numWindowsSeen;
                                }
                                windowStatus = true;
                            }
                        }

                        if (numWindowsSeen < minMatchingWindows) {
                            // a unique sequence
                            
                            itr.start(hpc);
                            while (itr.hasNext()) {
                                itr.next();
                                bf.add(hVals);
                            }

                            writer.write("seed" + Integer.toString(++id), seq);
                        }
                    }
                }
            }
        }
        
        writer.close();
        
        System.out.println("Bloom filter FPR:\t" + convertToRoundedPercent(bf.getFPR()) + " %");
        bf.destroy();
        
        System.out.println("before: " + NumberFormat.getInstance().format(numSeq) + 
                            "\tafter: " + NumberFormat.getInstance().format(id) +
                            " (" + convertToRoundedPercent(id/(float)numSeq) + " %)");
    }
    
    public static void main(String[] args) {
        //debug
    }
}
