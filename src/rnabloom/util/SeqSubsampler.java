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
import java.util.BitSet;
import java.util.HashSet;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.CountingBloomFilter;
import rnabloom.bloom.hash.CanonicalHashFunction;
import rnabloom.bloom.hash.CanonicalNTHashIterator;
import rnabloom.bloom.hash.HashFunction;
import rnabloom.bloom.hash.MinimizerHashIterator;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.io.FastaWriter;
import static rnabloom.util.Common.convertToRoundedPercent;
import static rnabloom.util.SeqUtils.compressHomoPolymers;

/**
 *
 * @author Ka Ming Nip
 */
public class SeqSubsampler {
    
    public static void minimizerBased(ArrayList<? extends BitSequence> seqs, String outFasta,
            long bfSize, int k, int w, int numHash, boolean stranded, boolean useHpcKmers,
            int maxNonMatchingChainLength, float minMatchingProportion,
            int maxMultiplicity) throws IOException {        

        int numSeq = seqs.size();
        
        HashFunction h = stranded ? new HashFunction(k) : new CanonicalHashFunction(k);
        
        MinimizerHashIterator itr = new MinimizerHashIterator(k, w, h.getHashIterator(1));
        CountingBloomFilter bf = new CountingBloomFilter(bfSize, numHash, h);
        
        FastaWriter fw  = new FastaWriter(outFasta, false);
        int seqID = 0;
        long[] hVals = new long[numHash];
        for (BitSequence s : seqs) {
            String seq = s.toString();
            String hpc = useHpcKmers ? compressHomoPolymers(seq) : seq;
            
            if (itr.start(hpc)) {
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
                
                if (maxConsecutiveMissing > maxNonMatchingChainLength || numMinimizersSeen < minMatchingProportion * numMinimizers) {
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
        
    public static void kmerBased(ArrayList<? extends BitSequence> seqs,
            String outSubsampleFasta, String outAllFasta,
            long bfSize, int k, int numHash, boolean stranded, 
            int maxMultiplicity, int maxEdgeClip, boolean verbose) throws IOException, InterruptedException {
                
        int numSeq = seqs.size();
        BitSet writeBits = new BitSet(numSeq);
        
        final int shift = k + 1;
        final int shiftD = k;
        final int shiftI = k + 2;
        final int missingChainThreshold = k + shift;
        int numSubsample = 0;
        float fpr = 0;

        if (stranded) {
            HashFunction h = new HashFunction(k);
            NTHashIterator hashItr = h.getHashIterator(1);
            CountingBloomFilter cbf = new CountingBloomFilter(bfSize, numHash, h);
            
            for (int seqIndex=0; seqIndex<numSeq; ++seqIndex) {
                String seq = seqs.get(seqIndex).toString();

                if (hashItr.start(seq)) {
                    int seqLen = seq.length();
                    boolean tooShort = seqLen < 3 * maxEdgeClip;
                    int missingChainLen = 0;
                    boolean write = false;
                    int numKmers = seqLen - k + 1;

                    // get hash value of all kmers
                    long[] kmerHashVals = new long[numKmers];
                    for (int i=0; i<numKmers; ++i) {
                        hashItr.next();
                        kmerHashVals[i] = hashItr.hVals[0];
                    }

                    int start = maxEdgeClip;
                    int end = numKmers - maxEdgeClip - shift;
                    if (tooShort) {
                        start = 0;
                        end = numKmers - shift;
                    }

                    // check whether this sequence has been seen enough by looking up multiplicities of k-mer pairs
                    for (int i=start; i<end; ++i) {
                        long pair = HashFunction.combineHashValues(kmerHashVals[i], kmerHashVals[i + shift]);

                        if (cbf.getCount(pair) <= maxMultiplicity) {
                            if (++missingChainLen >= missingChainThreshold) {
                                write = true;
                                break;
                            }
                        }
                        else {
                            missingChainLen = 0;
                        }
                    }

                    if (write) {
                        ++numSubsample;
                        writeBits.set(seqIndex);

                        // store k-mer pairs along this sequence
                        HashSet<Long> hashVals = new HashSet();

                        // k-mer pair gap size: 1
                        for (int i=start; i<end; ++i) {
                            hashVals.add(HashFunction.combineHashValues(kmerHashVals[i], kmerHashVals[i + shift]));
                        }

                        // k-mer pair gap size: 2 
                        for (int i=start; i<end - 1; ++i) {
                            hashVals.add(HashFunction.combineHashValues(kmerHashVals[i], kmerHashVals[i + shiftI]));
                        }

                        // k-mer pair gap size: 0
                        for (int i=start; i<end + 1; ++i) {
                            hashVals.add(HashFunction.combineHashValues(kmerHashVals[i], kmerHashVals[i + shiftD]));
                        }

                        hashVals.parallelStream().forEach(e -> {
                            cbf.increment(e);
                        });
                    }
                }
            }
            
            fpr = cbf.getFPR();
            cbf.destroy();
        }
        else {
            CanonicalHashFunction h = new CanonicalHashFunction(k);
            CanonicalNTHashIterator hashItr = (CanonicalNTHashIterator) h.getHashIterator(1);
            CountingBloomFilter cbf = new CountingBloomFilter(bfSize, numHash, h);
            
            for (int seqIndex=0; seqIndex<numSeq; ++seqIndex) {
                String seq = seqs.get(seqIndex).toString();

                if (hashItr.start(seq)) {
                    int seqLen = seq.length();
                    boolean tooShort = seqLen < 3 * maxEdgeClip;
                    int missingChainLen = 0;
                    boolean write = false;
                    int numKmers = seqLen - k + 1;

                    // get hash value of all kmers
                    long[] forwardKmerHashVals = new long[numKmers];
                    long[] reverseKmerHashVals = new long[numKmers];
                    for (int i=0; i<numKmers; ++i) {
                        hashItr.next();
                        forwardKmerHashVals[i] = hashItr.frhval[0];
                        reverseKmerHashVals[i] = hashItr.frhval[1];
                    }

                    int start = maxEdgeClip;
                    int end = numKmers - maxEdgeClip - shift;
                    if (tooShort) {
                        start = 0;
                        end = numKmers - shift;
                    }

                    // check whether this sequence has been seen enough by looking up multiplicities of k-mer pairs
                    for (int i=start; i<end; ++i) {
                        long pairF = HashFunction.combineHashValues(forwardKmerHashVals[i], forwardKmerHashVals[i + shift]);
                        long pairR = HashFunction.combineHashValues(reverseKmerHashVals[i + shift], reverseKmerHashVals[i]);
                        
                        if (cbf.getCount(Math.min(pairF, pairR)) <= maxMultiplicity) {
                            if (++missingChainLen >= missingChainThreshold) {
                                write = true;
                                break;
                            }
                        }
                        else {
                            missingChainLen = 0;
                        }
                    }

                    if (write) {
                        ++numSubsample;
                        writeBits.set(seqIndex);

                        // store k-mer pairs along this sequence
                        HashSet<Long> hashVals = new HashSet();

                        // k-mer pair gap size: 1
                        for (int i=start; i<end; ++i) {
                            long pairF = HashFunction.combineHashValues(forwardKmerHashVals[i], forwardKmerHashVals[i + shift]);
                            long pairR = HashFunction.combineHashValues(reverseKmerHashVals[i + shift], reverseKmerHashVals[i]);
                            hashVals.add(Math.min(pairF, pairR));
                        }

                        // k-mer pair gap size: 2 
                        for (int i=start; i<end - 1; ++i) {
                            long pairF = HashFunction.combineHashValues(forwardKmerHashVals[i], forwardKmerHashVals[i + shiftI]);
                            long pairR = HashFunction.combineHashValues(reverseKmerHashVals[i + shiftI], reverseKmerHashVals[i]);
                            hashVals.add(Math.min(pairF, pairR));
                        }

                        // k-mer pair gap size: 0
                        for (int i=start; i<end + 1; ++i) {
                            long pairF = HashFunction.combineHashValues(forwardKmerHashVals[i], forwardKmerHashVals[i + shiftD]);
                            long pairR = HashFunction.combineHashValues(reverseKmerHashVals[i + shiftD], reverseKmerHashVals[i]);
                            hashVals.add(Math.min(pairF, pairR));
                        }

                        hashVals.parallelStream().forEach(e -> {
                            cbf.increment(e);
                        });
                    }
                }
            }
            
            fpr = cbf.getFPR();
            cbf.destroy();
        }
        
        if (verbose) {
            System.out.println("Bloom filter FPR:\t" + convertToRoundedPercent(fpr) + " %");
            System.out.println("before: " + NumberFormat.getInstance().format(numSeq) + 
                                "\tafter: " + NumberFormat.getInstance().format(numSubsample) +
                                " (" + convertToRoundedPercent(numSubsample/(float)numSeq) + " %)");
        }
        
        FastaWriter subWriter  = new FastaWriter(outSubsampleFasta, false);
        FastaWriter allWriter  = new FastaWriter(outAllFasta, false);
        for (int seqIndex=0; seqIndex<numSeq; ++seqIndex) {
            String seq = seqs.get(seqIndex).toString();
            allWriter.write("a" + seqIndex, seq);
            if (writeBits.get(seqIndex)) {
                subWriter.write("s" + seqIndex, seq);
            }
        }
        subWriter.close();
        allWriter.close();       
    }
    
    public static void minimalSet(ArrayList<? extends BitSequence> seqs, String outFasta,
            long bfSize, int k, int numHash, boolean stranded, boolean useHpcKmers,
            int windowSize, int minMatchingWindows, float minMatchingProportion,
            int minSeqLen) throws IOException, InterruptedException {
        int numSeq = seqs.size();
        
        FastaWriter writer = new FastaWriter(outFasta, false);
        
        HashFunction h = stranded ? new HashFunction(k) : new CanonicalHashFunction(k);
        
        BloomFilter bf = new BloomFilter(bfSize, numHash, h);
        NTHashIterator itr = h.getHashIterator(numHash, k);
        long[] hVals = itr.hVals;
        int id = 0;
        final int minKmersNeededPerWindow = 2;
        
        for (BitSequence s : seqs) {
            if (s != null && s.length >= minSeqLen) {
                String seq = s.toString();
//                ArrayList<String> segments = trimLowComplexityRegions(s.toString(), 100);
//                for (String seq : segments) {
                    if (seq.length() >= minSeqLen) {
                        String hpc = useHpcKmers ? compressHomoPolymers(seq) : seq;

                        if (itr.start(hpc)) {
                            int numKmers = hpc.length() - k + 1;
                            int numWindows = numKmers/windowSize;
                            if (numKmers % windowSize > 0) {
                                ++numWindows;
                            }

//                            if (numWindows >= minMatchingWindows) {
                                int numWindowsSeen = 0;
                                int windowIndex = 0;
                                int numKmersSeenInWindow = 0;
                                while (itr.hasNext()) {
                                    itr.next();

                                    if (itr.getPos()/windowSize > windowIndex) {
                                        ++windowIndex;
                                        numKmersSeenInWindow = 0;
                                    }

                                    if (bf.lookup(hVals)) {
                                        if (++numKmersSeenInWindow == minKmersNeededPerWindow) {
                                            ++numWindowsSeen;
                                        }
                                    }
                                }

                                if (numWindowsSeen < Math.max(minMatchingWindows, Math.round(minMatchingProportion * numWindows))) {
                                //if (numWindowsSeen < minMatchingProportion * numWindows) {
                                    // a unique sequence

                                    itr.start(hpc);
                                    while (itr.hasNext()) {
                                        itr.next();
                                        bf.add(hVals);
                                    }

                                    writer.write("s" + Integer.toString(++id), seq);
                                }
//                            }
                        }
                    }
//                }
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
