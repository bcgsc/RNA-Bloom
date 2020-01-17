/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
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
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.NoSuchElementException;
import java.util.TreeSet;
import rnabloom.RNABloom.ReadPair;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.Kmer;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import static rnabloom.util.KmerBitsUtils.bitsToSeq;
import static rnabloom.util.KmerBitsUtils.seqToBits;
import static rnabloom.util.SeqUtils.*;

/**
 *
 * @author Ka Ming Nip
 */
public final class GraphUtils {

    public static boolean isValidSequence(ArrayList<Kmer> kmers, int k) {
        byte[] last = kmers.get(0).bytes;
        
        int numKmers = kmers.size();
        for (int i=1; i<numKmers; ++i) {
            byte[] current = kmers.get(i).bytes;
            
            if (!compareByteArrays(last, 1, k, current, 0, k-1)) {
                return false;
            }
            
            last = current;
        }
        
        return true;
    }
    
    private static boolean compareByteArrays(byte[] arr1, int start1, int end1, byte[] arr2, int start2, int end2) {
        int j = start2;
        for (int i=start1; i<end1; ++i) {
            if(arr1[i] != arr2[j++]) {
                return false;
            }
        }
        
        return true;
    }
    
    public static void printPairedKmersPositions(final ArrayList<Kmer> kmers, final BloomFilterDeBruijnGraph graph) {
        int d = graph.getFragPairedKmerDistance();
        int maxIndex = kmers.size() - 1 - d;
        
        System.out.println("FPK START");
        for (int i=0; i<=maxIndex; ++i) {
            int partnerIndex = i+d;
            if (graph.lookupKmerPair(kmers.get(i), kmers.get(partnerIndex))) {
                System.out.println(i + " " + partnerIndex + " " + getMinimumKmerCoverage(kmers, i, partnerIndex));
            }
        }
        System.out.println("FPK END");

        d = graph.getReadPairedKmerDistance();
        maxIndex = kmers.size() - 1 - d;
        
        System.out.println("RPK START");
        for (int i=0; i<=maxIndex; ++i) {
            int partnerIndex = i+d;
            if (graph.lookupReadKmerPair(kmers.get(i), kmers.get(partnerIndex))) {
                System.out.println(i + " " + partnerIndex + " " + getMinimumKmerCoverage(kmers, i, partnerIndex));
            }
        }        
        System.out.println("RPK END");
        
        System.out.println("R-SEG START");
        ArrayDeque<int[]> readRanges = breakWithReadPairedKmers(kmers, graph, 3);
        for (int[] r : readRanges) {
            System.out.println(r[0] + " " + r[1]);
        }
        System.out.println("R-SEG END");
        
        System.out.println("F-SEG START");
        ArrayDeque<int[]> fragRanges = breakWithFragPairedKmers(kmers, graph, 3);
        for (int[] r : fragRanges) {
            System.out.println(r[0] + " " + r[1]);
        }
        System.out.println("F-SEG END");
    }
    
    public static float getMinimumKmerCoverage(final Iterable<Kmer> kmers) {
        float min = Float.MAX_VALUE;
        
        float c;
        for (Kmer kmer : kmers) {
            c = kmer.count;
            if (c < min) {
                min = c;
            }
        }
        
        return min;
    }
    
    public static float getMinimumKmerCoverage(final ArrayList<Kmer> kmers, int start, int end) {
        float min = kmers.get(start).count;
        
        float c;
        for (int i=start+1; i<end; ++i) {
            c = kmers.get(i).count;
            if (c < min) {
                min = c;
            }
        }
        
        return min;
    }
    
    public static boolean areKmerCoverageAboveThreshold(final ArrayList<Kmer> kmers, int start, int end, float threshold) {
        float c;
        for (int i=start; i<end; ++i) {
            c = kmers.get(i).count;
            if (c < threshold) {
                return false;
            }
        }
        
        return true;
    }
    
    public static int getNumKmersAboveCoverageThreshold(final ArrayList<Kmer> kmers, int start, int end, float threshold) {
        float c;
        int num = 0;
        for (int i=start; i<end; ++i) {
            c = kmers.get(i).count;
            if (c >= threshold) {
                ++num;
            }
        }
        
        return num;
    }
    
    private static int getMinimumKmerCoverageIndexL2R(final ArrayList<Kmer> kmers, int start, int end) {
        // search from left to right
        
        float min = kmers.get(start).count;
        int index = start;
        
        float c;
        for (int i=start+1; i<end; ++i) {
            c = kmers.get(i).count;
            if (c < min) {
                min = c;
                index = i;
            }
        }
        
        return index;
    }
    
    private static int getMinimumKmerCoverageIndexR2L(final ArrayList<Kmer> kmers, int start, int end) {
        // search from right to left

        int index = end-1;
        float min = kmers.get(index).count;
        
        float c;
        for (int i=index-1; i>=start; --i) {
            c = kmers.get(i).count;
            if (c < min) {
                min = c;
                index = i;
            }
        }
        
        return index;
    }
    
    public static float getMedianKmerCoverage(final ArrayList<Kmer> kmers, int start, int end) {
        int range = end-start;
        
        float[] covs = new float[range];
        for (int i=0; i<range; ++i) {
            covs[i] = kmers.get(start+i).count;
        }

        return getMedian(covs);
    }
    
    public static float[] getMinMedMaxKmerCoverage(final ArrayList<Kmer> kmers) {        
        int len = kmers.size();
        float[] covs = new float[len];
        for (int i=0; i<len; ++i) {
            covs[i] = kmers.get(i).count;
        }

        return getMinMedMax(covs);
    }    
    
    public static float getMedianKmerCoverage(final Collection<Kmer> kmers) {
        int numKmers = kmers.size();
        
        Iterator<Kmer> itr = kmers.iterator();
        
        float[] counts = new float[numKmers];
        for (int i=0; i<numKmers; ++i) {
            counts[i] = itr.next().count;
        }
        Arrays.sort(counts);
        
        int halfNumKmers = numKmers/2;
        
        if (numKmers % 2 == 0) {
            return (counts[halfNumKmers] + counts[halfNumKmers-1])/2.0f;
        }
        
        return counts[halfNumKmers];
    }
    
    public static float getMedianKmerCoverage(final ArrayDeque<Kmer>... kmers) {
        int numKmers = 0;
        for (ArrayDeque<Kmer> c : kmers) {
            numKmers += c.size();
        }
        
        float[] counts = new float[numKmers];
        int i = -1;
        for (ArrayDeque<Kmer> c : kmers) {
            Iterator<Kmer> itr = c.iterator();
            while (itr.hasNext()) {
                counts[++i] = itr.next().count;
            }
        }
        
        Arrays.sort(counts);
        
        int halfNumKmers = numKmers/2;
        
        if (numKmers % 2 == 0) {
            return (counts[halfNumKmers] + counts[halfNumKmers-1])/2.0f;
        }
        
        return counts[halfNumKmers];
    }
    
    public static float getMaxMedianCoverageRight(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead) {
        ArrayDeque<Kmer> neighbors = source.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph);
        
        if (neighbors.isEmpty()) {
            if (lookahead > 0) {
                return 0;
            }
            else {
                return source.count;
            }
        }
        else {
            ArrayDeque<Kmer> path = new ArrayDeque<>(lookahead); 
            path.add(source);
            
            Kmer cursor = neighbors.removeFirst();
            path.add(cursor);

            ArrayDeque<ArrayDeque<Kmer>> frontier = new ArrayDeque<>(lookahead);
            frontier.add(neighbors);

            float bestCov = 0;

            while (!frontier.isEmpty()) {
                if (path.size() < lookahead) {
                    neighbors = cursor.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph);;
                    if (!neighbors.isEmpty()) {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        frontier.add(neighbors);
                        continue;
                    }
                }

                if (path.size() == lookahead) {
                    // we only calculate coverage if path is long enough
                    float pathCov = getMinimumKmerCoverage(path);
                    if (bestCov < pathCov) {
                        bestCov = pathCov;
                    }
                }

                while (!frontier.isEmpty()) {
                    neighbors = frontier.getLast();
                    path.removeLast();
                    if (neighbors.isEmpty()) {
                        frontier.removeLast();
                    }
                    else {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        break;
                    }
                }
            }
            
            return bestCov;
        }
    }
    
    public static float getMaxMedianCoverageRight(final BloomFilterDeBruijnGraph graph,
                                                final Kmer source,
                                                final int lookahead,
                                                final BloomFilter bf) {
        ArrayDeque<Kmer> neighbors = source.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph, bf);
        
        if (neighbors.isEmpty()) {
            if (lookahead > 0) {
                return 0;
            }
            else {
                return source.count;
            }
        }
        else {
            ArrayDeque<Kmer> path = new ArrayDeque<>(lookahead); 
            path.add(source);
            
            Kmer cursor = neighbors.removeFirst();
            path.add(cursor);

            ArrayDeque<ArrayDeque<Kmer>> frontier = new ArrayDeque<>(lookahead);
            frontier.add(neighbors);

            float bestCov = 0;

            while (!frontier.isEmpty()) {
                if (path.size() < lookahead) {
                    neighbors = cursor.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph, bf);
                    if (!neighbors.isEmpty()) {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        frontier.add(neighbors);
                        continue;
                    }
                }

                if (path.size() == lookahead) {
                    // we only calculate coverage if path is long enough
                    float pathCov = getMinimumKmerCoverage(path);
                    if (bestCov < pathCov) {
                        bestCov = pathCov;
                    }
                }

                while (!frontier.isEmpty()) {
                    neighbors = frontier.getLast();
                    path.removeLast();
                    if (neighbors.isEmpty()) {
                        frontier.removeLast();
                    }
                    else {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        break;
                    }
                }
            }
            
            return bestCov;
        }
    }

    public static float getMaxMedianCoverageLeft(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead) {
        ArrayDeque<Kmer> neighbors = source.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph);
        
        if (neighbors.isEmpty()) {
            if (lookahead > 0) {
                return 0;
            }
            else {
                return source.count;
            }
        }
        else {
            ArrayDeque<Kmer> path = new ArrayDeque<>(lookahead); 
            path.add(source);
            
            Kmer cursor = neighbors.removeFirst();
            path.add(cursor);

            ArrayDeque<ArrayDeque<Kmer>> frontier = new ArrayDeque<>(lookahead);
            frontier.add(neighbors);

            float bestCov = 0;

            while (!frontier.isEmpty()) {
                if (path.size() < lookahead) {
                    neighbors = cursor.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph);
                    if (!neighbors.isEmpty()) {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        frontier.add(neighbors);
                        continue;
                    }
                }

                if (path.size() == lookahead) {
                    // we only calculate coverage if path is long enough
                    float pathCov = getMinimumKmerCoverage(path);
                    if (bestCov < pathCov) {
                        bestCov = pathCov;
                    }
                }

                while (!frontier.isEmpty()) {
                    neighbors = frontier.getLast();
                    path.removeLast();
                    if (neighbors.isEmpty()) {
                        frontier.removeLast();
                    }
                    else {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        break;
                    }
                }
            }
            
            return bestCov;
        }
    }
    
    public static float getMaxMedianCoverageLeft(final BloomFilterDeBruijnGraph graph, 
                                                final Kmer source, 
                                                final int lookahead, 
                                                final BloomFilter bf) {
        ArrayDeque<Kmer> neighbors = source.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph, bf);
        
        if (neighbors.isEmpty()) {
            if (lookahead > 0) {
                return 0;
            }
            else {
                return source.count;
            }
        }
        else {
            ArrayDeque<Kmer> path = new ArrayDeque<>(lookahead); 
            path.add(source);
            
            Kmer cursor = neighbors.removeFirst();
            path.add(cursor);

            ArrayDeque<ArrayDeque<Kmer>> frontier = new ArrayDeque<>(lookahead);
            frontier.add(neighbors);

            float bestCov = 0;

            while (!frontier.isEmpty()) {
                if (path.size() < lookahead) {
                    neighbors = cursor.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph, bf);
                    if (!neighbors.isEmpty()) {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        frontier.add(neighbors);
                        continue;
                    }
                }

                if (path.size() == lookahead) {
                    // we only calculate coverage if path is long enough
                    float pathCov = getMinimumKmerCoverage(path);
                    if (bestCov < pathCov) {
                        bestCov = pathCov;
                    }
                }

                while (!frontier.isEmpty()) {
                    neighbors = frontier.getLast();
                    path.removeLast();
                    if (neighbors.isEmpty()) {
                        frontier.removeLast();
                    }
                    else {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        break;
                    }
                }
            }
            
            return bestCov;
        }
    }
    
    public static Kmer greedyExtendRightOnce(final BloomFilterDeBruijnGraph graph, final ArrayDeque<Kmer> candidates, final int lookahead) {
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.peek();
            }
            else {
                float bestCov = -1;
                Kmer bestKmer = null;
                for (Kmer kmer : candidates) {
                    float c = getMaxMedianCoverageRight(graph, kmer, lookahead);
                    if (c > bestCov) {
                        bestKmer = kmer;
                        bestCov = c;
                    }
                    else if (c == bestCov && kmer.count > bestKmer.count) {
                        bestKmer = kmer;
                    }
                }
                return bestKmer;
            }
        }
    }
    
    public static Kmer greedyExtendRightOnce(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead) {
        return greedyExtendRightOnce(graph, source.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph), lookahead);
    }
    
    public static Kmer greedyExtendRightOnce(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead, BloomFilter bf) {
        return greedyExtendRightOnce(graph, source.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph, bf), lookahead, bf);
    }
    
    public static Kmer greedyExtendRightOnce(final BloomFilterDeBruijnGraph graph,
                                            final ArrayDeque<Kmer> candidates,
                                            final int lookahead,
                                            final BloomFilter bf) {
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.peek();
            }
            else {
                float bestCov = -1;
                Kmer bestKmer = null;
                for (Kmer kmer : candidates) {
                    float c = getMaxMedianCoverageRight(graph, kmer, lookahead, bf);
                    if (c > bestCov) {
                        bestKmer = kmer;
                        bestCov = c;
                    }
                    else if (c == bestCov && kmer.count > bestKmer.count) {
                        bestKmer = kmer;
                    }
                }
                return bestKmer;
            }
        }
    }
    
    public static Kmer greedyExtendLeftOnce(final BloomFilterDeBruijnGraph graph, final ArrayDeque<Kmer> candidates, final int lookahead) {
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.peek();
            }
            else {
                float bestCov = -1;
                Kmer bestKmer = null;
                for (Kmer kmer : candidates) {
                    float c = getMaxMedianCoverageLeft(graph, kmer, lookahead);
                    if (c > bestCov) {
                        bestKmer = kmer;
                        bestCov = c;
                    }
                    else if (c == bestCov && kmer.count > bestKmer.count) {
                        bestKmer = kmer;
                    }
                }
                return bestKmer;
            }
        }
    }
    
    public static Kmer greedyExtendLeftOnce(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead) {
        return greedyExtendLeftOnce(graph, source.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph), lookahead);
    }

    public static Kmer greedyExtendLeftOnce(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead, BloomFilter bf) {
        return greedyExtendLeftOnce(graph, source.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph, bf), lookahead, bf);
    }
    
    public static Kmer greedyExtendLeftOnce(final BloomFilterDeBruijnGraph graph, 
                                            final ArrayDeque<Kmer> candidates, 
                                            final int lookahead, 
                                            final BloomFilter bf) {
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.peek();
            }
            else {
                float bestCov = -1;
                Kmer bestKmer = null;
                for (Kmer kmer : candidates) {
                    float c = getMaxMedianCoverageLeft(graph, kmer, lookahead, bf);
                    if (c > bestCov) {
                        bestKmer = kmer;
                        bestCov = c;
                    }
                    else if (c == bestCov && kmer.count > bestKmer.count) {
                        bestKmer = kmer;
                    }
                }
                return bestKmer;
            }
        }
    }
    
    public static boolean containsAllKmers(final BloomFilter bf,
                                final ArrayList<Kmer> kmers) {
        if (kmers.isEmpty()) {
            return false;
        }
        
        for (Kmer kmer : kmers) {
            if (kmer == null || !bf.lookup(kmer.getHash())) {
                return false;
            }
        }
        
        return true;
    }

    private static class Sequence implements Comparable<Object> {
        int length;
        BitSet bits;
        
        public Sequence(String seq) {
            this.length = seq.length();
            this.bits = seqToBits(seq);
        }
        
        @Override
        public String toString() {
            return bitsToSeq(this.bits, this.length);
        }
        
        @Override
        public int compareTo(Object other) {
            return ((Sequence) other).length - this.length;
        }
    }
    
    public static int reduceRedundancy(final String inFasta,
                                        final String outFasta,
                                        final BloomFilterDeBruijnGraph graph,
                                        final BloomFilter bf,
                                        final int lookahead,
                                        final int maxIndelSize,
                                        final int maxTipLength,
                                        final float percentIdentity) throws IOException {
        ArrayList<Sequence> seqs = new ArrayList<>();
        
        // read entire FASTA and store all sequences
        FastaReader fr = new FastaReader(inFasta);
        try {
            while (true) {
                seqs.add(new Sequence(fr.next()));
            }
        }
        catch (NoSuchElementException e) {
            //end of file
        }
        fr.close();
        
        // sort sequences by length
        Collections.sort(seqs);

        // remove redundant sequences
        int cid = 0;
        FastaWriter fw = new FastaWriter(outFasta, false);
        for (Sequence s : seqs) {
            int len = s.length;
            
            String seq = s.toString();
            ArrayList<Kmer> kmers = graph.getKmers(seq);
            
            if (!represented(kmers, graph, bf, lookahead, maxIndelSize, maxTipLength, percentIdentity)) {
                // insert kmers into Bloom filter
                for (Kmer kmer : kmers) {
                    bf.add(kmer.getHash());
                }
                
                // write to file
                fw.write(++cid+" l="+len, seq);
            }
        }
        fw.close();
        
        return seqs.size() - cid;
    }
    
    public static boolean represented(final String[] kmers, final BloomFilter bf) {
        for (String kmer : kmers) {
            if (!bf.lookup(kmer)) {
                return false;
            }
        }
        
        return true;
    }
    
    public static boolean represented(final ArrayList<Kmer> kmers,
                                    final BloomFilterDeBruijnGraph graph,
                                    final BloomFilter bf,
                                    final int lookahead,
                                    final int maxIndelSize,
                                    final int maxEdgeClipLength,
                                    final float percentIdentity) {
        int numKmers = kmers.size();
        int maxIndex = numKmers - 1;
        
        int k = graph.getK();
        int maxNumBubbleKmers = graph.getReadPairedKmerDistance()+k;
        
        int lastRepresentedKmerFoundIndex = -1;
        
        for (int i=0; i<numKmers; ++i) {
            
            if (bf.lookup(kmers.get(i).getHash())) {
                int startIndex = i;
                int endIndex = i;
                for (int j=i+1; j<numKmers; ++j) {
                    if (bf.lookup(kmers.get(j).getHash())) {
                        endIndex = j;
                    }
                    else {
                        break;
                    }
                }
                
                int assembledRange = endIndex - startIndex + 1;
                
                if (assembledRange >= lookahead) {                    
                    if (startIndex > 0) {
                                                
                        if (lastRepresentedKmerFoundIndex < 0) {                            
                            if (startIndex >= maxEdgeClipLength || hasDepthLeft(kmers.get(0), graph, maxEdgeClipLength-startIndex)) {
                                // check left edge kmers
                                ArrayDeque<Kmer> testEdgeKmers = greedyExtendLeft(graph, kmers.get(startIndex), lookahead, startIndex, bf);
                                if (testEdgeKmers.size() != startIndex ||
                                        getPercentIdentity(graph.assemble(testEdgeKmers), graph.assemble(kmers, 0, startIndex)) < percentIdentity) {
                                    return false;
                                }
                            }
                        }
                        else {
                            // check gap kmers
                            int expectedPathLen = startIndex-lastRepresentedKmerFoundIndex-1;

                            if (expectedPathLen > maxNumBubbleKmers) {
                                return false;
                            }
                            
                            int numMissing = k-expectedPathLen;
                            
                            int left = lastRepresentedKmerFoundIndex;
                            int right = startIndex;
                            
                            if (numMissing > 0) {
                                for (int j=0; j<numMissing; ++j) {
                                    if (left == 0 || !bf.lookup(kmers.get(--left).getHash())) {
                                        break;
                                    }
                                }

                                for (int j=0; j<numMissing; ++j) {
                                    if (right == maxIndex || !bf.lookup(kmers.get(++right).getHash())) {
                                        break;
                                    }
                                }
                                
                                expectedPathLen = right - left - 1;
                            }
                            
                            ArrayDeque<Kmer> testPathKmers = getMaxCoveragePath(graph, kmers.get(left), kmers.get(right), expectedPathLen+maxIndelSize, lookahead, bf);
                            if (testPathKmers == null) {
                                return false;
                            }

                            int testPathLen = testPathKmers.size();

                            if ((testPathLen < expectedPathLen-maxIndelSize ||
                                    testPathLen > expectedPathLen+maxIndelSize ||
                                        getPercentIdentity(graph.assemble(testPathKmers), graph.assemble(kmers, left+1, right)) < percentIdentity)) {
                                return false;
                            }
                        }
                    }
                        
                    lastRepresentedKmerFoundIndex = endIndex;
                }
                
                i = endIndex;
            }
        }
        
        if (lastRepresentedKmerFoundIndex >= 0) {
            if (lastRepresentedKmerFoundIndex < maxIndex) {
                // check right edge kmers
                int expectedLen = numKmers-lastRepresentedKmerFoundIndex-1;
                if (expectedLen >= maxEdgeClipLength || hasDepthRight(kmers.get(maxIndex), graph, maxEdgeClipLength-expectedLen)) {                    
                    ArrayDeque<Kmer> testEdgeKmers = greedyExtendRight(graph, kmers.get(lastRepresentedKmerFoundIndex), lookahead, expectedLen, bf);
                    if (testEdgeKmers.size() != expectedLen ||
                            getPercentIdentity(graph.assemble(testEdgeKmers), graph.assemble(kmers, lastRepresentedKmerFoundIndex+1, numKmers)) < percentIdentity) {
                        return false;
                    }
                }
            }
        }
        else {
            return false;
        }
        
        return true;
    }
    
    public static boolean hasValidPath(BloomFilterDeBruijnGraph graph,
                                       Kmer left, 
                                       Kmer right, 
                                       BloomFilter bf, 
                                       int lowerBound, 
                                       int upperBound) {
        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<Kmer> frontier = new ArrayDeque<>();
        frontier.addAll(left.getSuccessors(k, numHash, graph, bf));
        
        HashSet<Kmer> kmersInFrontier = new HashSet<>();
        ArrayDeque<Kmer> newFrontier;
        for (int i=1; i<lowerBound; ++i) {
            kmersInFrontier.clear();
            newFrontier = new ArrayDeque<>();
            for (Kmer kmer : frontier) {
                for (Kmer s : kmer.getSuccessors(k, numHash, graph, bf)) {
                    if (!kmersInFrontier.contains(s)) { 
                        newFrontier.add(s);
                        kmersInFrontier.add(s);
                    }
                }
            }
            
            if (newFrontier.isEmpty()) {
                return false;
            }
            
            frontier = newFrontier;
        }
        
        for (int i=lowerBound; i<=upperBound; ++i) {
            kmersInFrontier.clear();
            newFrontier = new ArrayDeque<>();
            for (Kmer kmer : frontier) {
                if (kmer.equals(right)) {
                    return true;
                }
                newFrontier.add(kmer);
                for (Kmer s : kmer.getSuccessors(k, numHash, graph, bf)) {
                    if (!kmersInFrontier.contains(s)) { 
                        newFrontier.add(s);
                        kmersInFrontier.add(s);
                    }
                }
            }
            
            if (newFrontier.isEmpty()) {
                return false;
            }
            
            frontier = newFrontier;
        }
        
        for (Kmer kmer : frontier) {
            if (kmer.equals(right)) {
                return true;
            }
        }
        
        return false;
    }
    
    public static ArrayList<Kmer> join(BloomFilterDeBruijnGraph graph, 
                                        ArrayList<Kmer> leftKmers, 
                                        ArrayList<Kmer> rightKmers, 
                                        int bound, 
                                        int lookahead, 
                                        float maxCovGradient,
                                        int maxTipLen,
                                        int maxIndelLen,
                                        float minPercentIdentity,
                                        float minKmerCov) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        //int pkd = graph.getReadPairedKmerDistance();
        
        HashSet<Kmer> leftKmersSet = new HashSet<>(leftKmers);
        HashSet<Kmer> rightKmersSet = new HashSet<>(rightKmers);
        

        Kmer right = rightKmers.get(0);
        
        float leftCoverageThreshold = Math.max(minKmerCov, getMinimumKmerCoverage(leftKmers) * maxCovGradient);
        
        ArrayList<Kmer> leftPath = new ArrayList<>(leftKmers);
        
        int maxSize = bound + leftPath.size();
        while (leftPath.size() < maxSize) {            
            ArrayDeque<Kmer> neighbors = leftPath.get(leftPath.size()-1).getSuccessors(k, numHash, graph, minKmerCov);

            // filter by coverage
            if (neighbors.isEmpty()) {
                break;
            }
            else {
                Kmer best = null;

                if (neighbors.size() == 1) {
                    best = neighbors.pop();
                    leftCoverageThreshold = Math.max(minKmerCov, Math.min(leftCoverageThreshold, best.count*maxCovGradient));
                }
                else {
                    Iterator<Kmer> itr = neighbors.iterator();
                    while (itr.hasNext()) {
                        if (itr.next().count < leftCoverageThreshold) {
                            itr.remove();
                        }
                    }

                    if (neighbors.size() == 1) {
                        best = neighbors.pop();
                        leftCoverageThreshold = Math.max(minKmerCov, Math.min(leftCoverageThreshold, best.count*maxCovGradient));
                    }
                }

                if (best != null) {                    
                    if (rightKmersSet.contains(best)) {
                        if (best.equals(right)) {
                            leftPath.addAll(rightKmers);
                            return leftPath;
                        }
                        else {
                            int bestIndex = rightKmers.indexOf(best);

                            float pathCov = getMinimumKmerCoverage(leftPath, Math.max(0, leftPath.size()-bestIndex), leftPath.size());
                            float danglingCov = getMinimumKmerCoverage(rightKmers, 0, bestIndex);

                            if (danglingCov < pathCov) {
                                ArrayList<Kmer> fragmentKmers = new ArrayList<>(leftPath.size() + rightKmers.size() - bestIndex);
                                fragmentKmers.addAll(leftPath);
                                for (int i=bestIndex; i<rightKmers.size(); ++i) {
                                    fragmentKmers.add(rightKmers.get(i));
                                }

                                return fragmentKmers;
                            }
                        }
                    }

                    if (isHomopolymer(best.bytes) || leftKmersSet.contains(best)) {
                        break;
                    }
                    else {
                        leftPath.add(best);
                        leftKmersSet.add(best);
                    }
                }
                else {
                    ArrayDeque<Kmer> e = extendRightSE(leftPath, graph, maxTipLen, minKmerCov);

                    if (e == null || e.isEmpty()) {
                        break;
                    }

                    boolean atLoop = false;
                    
                    for (Kmer current : e) {
                        if (rightKmersSet.contains(current)) {
                            if (current.equals(right)) {
                                leftPath.addAll(rightKmers);
                                return leftPath;
                            }
                            else {
                                int bestIndex = rightKmers.indexOf(current);

                                float pathCov = getMinimumKmerCoverage(leftPath, Math.max(0, leftPath.size()-bestIndex), leftPath.size());
                                float danglingCov = getMinimumKmerCoverage(rightKmers, 0, bestIndex);

                                if (danglingCov < pathCov) {
                                    ArrayList<Kmer> fragmentKmers = new ArrayList<>(leftPath.size() + rightKmers.size() - bestIndex);
                                    fragmentKmers.addAll(leftPath);
                                    for (int i=bestIndex; i<rightKmers.size(); ++i) {
                                        fragmentKmers.add(rightKmers.get(i));
                                    }

                                    return fragmentKmers;
                                }
                            }
                        }
                        
                        if (isHomopolymer(current.bytes)) {
                            atLoop = true;
                            break;
                        }
                        else {
                            leftPath.add(current);
                            leftKmersSet.add(current);

                            float c = maxCovGradient * current.count;
                            if (c < leftCoverageThreshold) {
                                leftCoverageThreshold = Math.max(c, minKmerCov);
                            }
                        }
                    }
                    
                    if (atLoop) {
                        break;
                    }
                }
            }
        }
        
        Kmer left = leftKmers.get(leftKmers.size()-1);
        float rightCoverageThreshold = Math.max(minKmerCov, getMinimumKmerCoverage(rightKmers) * maxCovGradient);
        
        ArrayList<Kmer> rightPath = new ArrayList<>(rightKmers);
        Collections.reverse(rightPath);
        
        maxSize = bound + rightPath.size();
        while (rightPath.size() < maxSize) {
            ArrayDeque<Kmer> neighbors = rightPath.get(rightPath.size()-1).getPredecessors(k, numHash, graph, minKmerCov);

            // filter by coverage
            if (neighbors.isEmpty()) {
                break;
            }
            else {
                Kmer best = null;

                if (neighbors.size() == 1) {
                    best = neighbors.pop();
                    rightCoverageThreshold = Math.max(minKmerCov, Math.min(rightCoverageThreshold, best.count*maxCovGradient));
                }
                else {
                    Iterator<Kmer> itr = neighbors.iterator();
                    while (itr.hasNext()) {
                        if (itr.next().count < rightCoverageThreshold) {
                            itr.remove();
                        }
                    }

                    if (neighbors.size() == 1) {
                        best = neighbors.pop();
                        rightCoverageThreshold = Math.max(minKmerCov, Math.min(rightCoverageThreshold, best.count*maxCovGradient));
                    }
                }

                if (best != null) {
                    if (leftKmersSet.contains(best)) {
                        if (best.equals(left)) {
                            ArrayList<Kmer> fragmentKmers = new ArrayList<>(leftKmers.size() + rightPath.size());
                            fragmentKmers.addAll(leftKmers);
                            for (int i=rightPath.size()-1; i>=0; --i) {
                                fragmentKmers.add(rightPath.get(i));
                            }

                            return fragmentKmers;
                        }
                        else {
                            for (int i=leftPath.size()-1; i>=0; --i) {
                                if (leftPath.get(i).equals(best)) {
                                    break;
                                }
                                else {
                                    leftPath.remove(i);
                                }
                            }

                            for (int i=rightPath.size()-1; i>=0; --i) {
                                leftPath.add(rightPath.get(i));
                            }

                            return leftPath;                         
                        }
                    }
                    else {
                        if (isHomopolymer(best.bytes) || rightKmersSet.contains(best)) {
                            return null;
                        }
                        else {
                            rightPath.add(best);
                            rightKmersSet.add(best);
                        }
                    }
                }
                else {
                    ArrayDeque<Kmer> e = extendLeftSE(rightPath, graph, maxTipLen, minKmerCov);

                    if (e == null || e.isEmpty()) {
                        break;
                    }

                    for (Kmer current : e) {
                        if (leftKmersSet.contains(current)) {
                            for (int i=leftPath.size()-1; i>=0; --i) {
                                if (leftPath.get(i).equals(current)) {
                                    break;
                                }
                                else {
                                    leftPath.remove(i);
                                }
                            }

                            for (int i=rightPath.size()-1; i>=0; --i) {
                                leftPath.add(rightPath.get(i));
                            }

                            return leftPath;
                        }
                        else if (isHomopolymer(current.bytes)) {
                            return null;
                        }
                        else {
                            rightPath.add(current);
                            rightKmersSet.add(current);

                            float c = maxCovGradient * current.count;
                            if (c < rightCoverageThreshold) {
                                rightCoverageThreshold = Math.max(c, minKmerCov);
                            }
                        }
                    }
                }
            }
        }
        
        return null;
    }
    
    public static ArrayList<Kmer> getSimilarCoveragePath(BloomFilterDeBruijnGraph graph, 
                                                        ArrayList<Kmer> leftKmers, 
                                                        ArrayList<Kmer> rightKmers, 
                                                        int bound, 
                                                        int lookahead, 
                                                        float maxCovGradient,
                                                        boolean rescue,
                                                        float minKmerCov) {
        
        HashSet<Kmer> leftKmersSet = new HashSet<>(leftKmers);
        HashSet<Kmer> rightKmersSet = new HashSet<>(rightKmers);
        
        Kmer left = leftKmers.get(leftKmers.size()-1);
        Kmer right = rightKmers.get(0);
        
        float leftCoverageThreshold = getMinimumKmerCoverage(leftKmers) * maxCovGradient;
        float rightCoverageThreshold = getMinimumKmerCoverage(rightKmers) * maxCovGradient;
        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        HashSet<Kmer> leftPathKmers = new HashSet<>(bound);
        
        /* extend right */
        ArrayDeque<Kmer> leftPath = new ArrayDeque<>(bound);
        Kmer best;
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        
        best = left;

        for (int depth=0; depth < bound; ++depth) {
            best.getSuccessors(k, numHash, graph, neighbors, minKmerCov);

            if (neighbors.isEmpty()) {
                break;
            }
            else {
                if (neighbors.size() == 1) {
                    best = neighbors.pop();
                }
                else {
                    Iterator<Kmer> itr = neighbors.iterator();
                    while (itr.hasNext()) {
                        if (itr.next().count < leftCoverageThreshold) {
                            itr.remove();
                        }
                    }
                    
                    if (neighbors.isEmpty()) {
                        break;
                    }
                    else if (neighbors.size() == 1) {
                        best = neighbors.pop();
                        
                        float c = maxCovGradient * best.count;
                        if (c < leftCoverageThreshold) {
                            leftCoverageThreshold = c;
                        }
                    }
                    else {                        
//                        System.out.println(">L\n" + graph.assemble(leftKmers));
//                        System.out.println(">L_path\n" + graph.assemble(leftPath));
//                        System.out.println(">R\n" + graph.assemble(rightKmers));
                        break;
                    }
                }

                if (rightKmersSet.contains(best)) {
                    if (best.equals(right)) {
                        ArrayList<Kmer> fragmentKmers = new ArrayList<>(leftKmers.size() + leftPath.size() + rightKmers.size());
                        fragmentKmers.addAll(leftKmers);
                        fragmentKmers.addAll(leftPath);
                        fragmentKmers.addAll(rightKmers);
                        
                        return fragmentKmers;
                    }
                    else {
                        int bestIndex = rightKmers.indexOf(best);
                        
                        float pathCov = -1;
                        Iterator<Kmer> itr = leftPath.descendingIterator();
                        for (int i=0; i<bestIndex; ++i) {
                            if (itr.hasNext()) {
                                float c = itr.next().count;
                                if (c < pathCov || pathCov == -1) {
                                    pathCov = c;
                                }
                            }
                            else {
                                break;
                            }
                        }
                        
                        float danglingCov = getMinimumKmerCoverage(rightKmers, 0, bestIndex);
                        
                        if (danglingCov < pathCov) {
                            ArrayList<Kmer> fragmentKmers = new ArrayList<>(leftKmers.size() + leftPath.size() + rightKmers.size() - bestIndex);
                            fragmentKmers.addAll(leftKmers);
                            fragmentKmers.addAll(leftPath);
                            for (int i=bestIndex; i<rightKmers.size(); ++i) {
                                fragmentKmers.add(rightKmers.get(i));
                            }

//                        System.out.println(">L\n" + graph.assemble(leftKmers));
//                        System.out.println(">R\n" + graph.assemble(rightKmers));
//                        System.out.println(">F\n" + graph.assemble(fragmentKmers));
                            
                            return fragmentKmers;
                        }
                        else {
                            break;
                        }
                    }
                }
                else {
                    if (leftPathKmers.contains(best) || leftKmersSet.contains(best)) {
                        break;
                    }
                    else {
                        leftPathKmers.add(best);
                        leftPath.add(best);
                    }
                }
            }
        }
        
        neighbors.clear();
        HashSet<Kmer> rightPathKmers = new HashSet<>(bound);
        
        /* not connected, search from right */
        ArrayDeque<Kmer> rightPath = new ArrayDeque<>(bound);
        best = right;
        for (int depth=0; depth < bound; ++depth) {
            best.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
            
            if (neighbors.isEmpty()) {
                break;
            }
            else {
                if (neighbors.size() == 1) {
                    best = neighbors.pop();
                }
                else {
                    Iterator<Kmer> itr = neighbors.iterator();
                    while (itr.hasNext()) {
                        if (itr.next().count < rightCoverageThreshold) {
                            itr.remove();
                        }
                    }
                    
                    if (neighbors.isEmpty()) {
                        break;
                    }
                    else if (neighbors.size() == 1) {
                        best = neighbors.pop();
                        
                        float c = maxCovGradient * best.count;
                        if (c < rightCoverageThreshold) {
                            rightCoverageThreshold = c;
                        }
                    }
                    else {
//                        System.out.println(">L\n" + graph.assemble(leftKmers));
//                        System.out.println(">L_path\n" + graph.assemble(leftPath));
//                        System.out.println(">R_path\n" + graph.assemble(rightPath));
//                        System.out.println(">R\n" + graph.assemble(rightKmers));
                        break;
                    }
                }
                
                if (leftKmersSet.contains(best)) {
                    if (best.equals(left)) {
                        ArrayList<Kmer> fragmentKmers = new ArrayList<>(leftKmers.size() + rightPath.size() + rightKmers.size());
                        fragmentKmers.addAll(leftKmers);
                        fragmentKmers.addAll(rightPath);
                        fragmentKmers.addAll(rightKmers);
                        
                        return fragmentKmers;
                    }
                    else {
//                        System.out.println(">L\n" + graph.assemble(leftKmers));
//                        System.out.println(">L_path\n" + graph.assemble(leftPath));
//                        System.out.println(">R_path\n" + graph.assemble(rightPath));
//                        System.out.println(">R\n" + graph.assemble(rightKmers));
                        
                        int bestIndex = leftKmers.indexOf(best);
                        
                        float danglingCov = getMinimumKmerCoverage(leftKmers, bestIndex, leftKmers.size());
                        
                        float pathCov = -1;
                        Iterator<Kmer> itr = rightPath.iterator();
                        for (int i=leftKmers.size(); i>bestIndex; --i) {
                            if (itr.hasNext()) {
                                float c = itr.next().count;
                                if (c < pathCov || pathCov == -1) {
                                    pathCov = c;
                                }
                            }
                            else {
                                break;
                            }
                        }
                        
                        if (danglingCov < pathCov) {
                            ArrayList<Kmer> fragmentKmers = new ArrayList<>(bestIndex + rightPath.size() + rightKmers.size());
                            for (int i=0; i<=bestIndex; ++i) {
                                fragmentKmers.add(leftKmers.get(i));
                            }                            
                            fragmentKmers.addAll(rightPath);
                            fragmentKmers.addAll(rightKmers);
                            
//                            System.out.println(">L\n" + graph.assemble(leftKmers));
//                            System.out.println(">R\n" + graph.assemble(rightKmers));
//                            System.out.println(">f\n" + graph.assemble(fragmentKmers));
                            
                            return fragmentKmers;
                        }
                        else {
                            return null;
                        }
                    }
                }
                else {
                    if (rightPathKmers.contains(best) || rightKmersSet.contains(best)) {
                        return null;
                    }
                    else if (leftPathKmers.contains(best)) {
                        /* right path intersects the left path */
                        
//                        if (graph.isLowComplexity(best) &&
//                                (best.hasAtLeastXSuccessors(k, numHash, graph, 2) ||
//                                    best.hasAtLeastXPredecessors(k, numHash, graph, 2))) {
//                            return null;
//                        }
                        
                        ArrayList<Kmer> fragmentKmers = new ArrayList<>(leftKmers.size() + rightPath.size() + rightKmers.size());
                        fragmentKmers.addAll(leftKmers);
                        
                        for (Kmer kmer : leftPath) {
                            fragmentKmers.add(kmer);
                            
                            if (kmer.equals(best)) {
                                break;
                            }
                        }
                        
                        fragmentKmers.addAll(rightPath);
                        fragmentKmers.addAll(rightKmers);

                        return fragmentKmers;
                    }
                    else {
                        rightPathKmers.add(best);
                        rightPath.addFirst(best);
                    }
                }
            }
        }
        
        /*
        if (rescue) {
            // left path does not intersect right path
            
            for (Kmer kmer : leftPath) {
                ArrayDeque<Kmer> successors = kmer.getSuccessors(k, numHash, graph);
                if (successors.size() > 1) {
                    for (Kmer s : successors) {
                        if (!leftPathKmers.contains(s) ) {
                            if (rightPathKmers.contains(s)) {
                                ArrayDeque<Kmer> path = new ArrayDeque<>(leftPath.size() + rightPath.size());

                                // fill with kmers in left path
                                for (Kmer tmp : leftPath) {
                                    path.add(tmp);
                                    if (tmp.equals(kmer)) {
                                        break;
                                    }
                                }
                                
                                path.add(s);
                                
                                // fill with kmers in right path
                                boolean found = false;
                                for (Kmer tmp : rightPath) {
                                    if (found) {
                                        path.add(tmp);
                                    }
                                    else if (tmp.equals(s)) {
                                        found = true;
                                    }
                                }

                                return path;
                            }
                            else {
                                // perform a bounded greedy extension (depth = k) to connect left path to right path
                                ArrayDeque<Kmer> extension = greedyExtendRight(graph, s, lookahead, k);
                                for (Kmer e : extension) {
//                                    if (graph.isLowComplexity(e)) {
//                                        break;
//                                    }
                                    
                                    if (rightPathKmers.contains(e)) {
                                        ArrayDeque<Kmer> path = new ArrayDeque<>(leftPath.size() + rightPath.size());

                                        // fill with kmers in left path
                                        for (Kmer tmp : leftPath) {
                                            path.add(tmp);
                                            if (tmp.equals(kmer)) {
                                                break;
                                            }
                                        }

                                        path.add(s);
                                        
                                        // fill with kmers in extension
                                        for (Kmer tmp : extension) {
                                            path.add(tmp);
                                            if (tmp.equals(e)) {
                                                break;
                                            }
                                        }
                                        
                                        // fill with kmers in right path
                                        boolean found = false;
                                        for (Kmer tmp : rightPath) {
                                            if (found) {
                                                path.add(tmp);
                                            }
                                            else if (tmp.equals(e)) {
                                                found = true;
                                            }
                                        }
                                        
                                        return path;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        */
        
        return null;
    }
    
    public static ArrayDeque<Kmer> findPath(BloomFilterDeBruijnGraph graph, Kmer left, Kmer right, int bound, int lookahead, float minKmerCov) {
        if (!graph.isLowComplexity(left) && !graph.isLowComplexity(right)) {
            int k = graph.getK();
            int numHash = graph.getMaxNumHash();

            ArrayDeque<Kmer> rightExtension = new ArrayDeque<>();

            ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
            for (int i=0; i<bound; ++i) {
                right.getPredecessors(k, numHash, graph, neighbors, minKmerCov);

                if (neighbors.size() == 1) {
                    Kmer kmer = neighbors.pop();

                    if (left.equals(kmer)) {
                        return rightExtension;
                    }

                    rightExtension.addFirst(kmer);
                }
                else {
                    break;
                }
            }

            if (!rightExtension.isEmpty()) {
                bound -= rightExtension.size();
                right = rightExtension.getFirst();
            }

            // data structure to store visited kmers at defined depth
            HashSet<Kmer> visitedBranchingKmers = new HashSet<>();

            int depth = 0;

            ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
            branchesStack.add(getSuccessorsRanked(left, graph, lookahead));

            ArrayDeque<Kmer> extension = new ArrayDeque<>();
            HashSet<Kmer> extensionKmers = new HashSet<>();

            while (!branchesStack.isEmpty()) {
                LinkedList<Kmer> branches = branchesStack.getLast();

                if (branches.isEmpty()) {
                    extensionKmers.remove(extension.pollLast());
                    branchesStack.removeLast();
                    --depth;
                }
                else {
                    Kmer cursor = branches.pop();

                    if (cursor.equals(right)) {
                        if (!rightExtension.isEmpty()) {
                            extension.addAll(rightExtension);
                        }

                        return extension;
                    }

                    if (depth < bound && !extensionKmers.contains(cursor) && !graph.isLowComplexity(cursor)) {
                        if (cursor.hasAtLeastXPredecessors(k, numHash, graph, 2)) {
                            // these kmers may be visited from an alternative branch upstream
                            if (visitedBranchingKmers.add(cursor)) {
                                branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursor);
                                ++depth;
                            }
                        }
                        else {
                            branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursor);
                            ++depth;
                        }
                    }
                }
            }
        }
        
        return null;
    }
    
    /**
     * 
     * @param graph
     * @param left
     * @param right
     * @param bound
     * @param lookahead
     * @param minKmerCov
     * @return 
     */
    public static ArrayDeque<Kmer> getMaxCoveragePath(BloomFilterDeBruijnGraph graph, Kmer left, Kmer right, int bound, int lookahead, float minKmerCov) {
               
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        HashSet<Kmer> leftPathKmers = new HashSet<>(bound);
        
        /* extend right */
        ArrayDeque<Kmer> leftPath = new ArrayDeque<>(bound);
        Kmer best;
        ArrayDeque<Kmer> neighbors;
        
        best = left;

        for (int depth=0; depth < bound; ++depth) {
            neighbors = best.getSuccessors(k, numHash, graph, minKmerCov);
//            best.successors = null; // clear cache

            if (neighbors.isEmpty()) {
                break;
            }
            else {
                if (neighbors.size() == 1) {
                    best = neighbors.peek();
                }
                else {
                    best = greedyExtendRightOnce(graph, neighbors, lookahead);
                }

                if (best.equals(right)) {
                    return leftPath;
                }
                else {
                    if (leftPathKmers.contains(best)) {
                        break;
                    }
                    else {
                        leftPathKmers.add(best);
                        leftPath.add(best);
                    }
                }
            }
        }
        
        HashSet<Kmer> rightPathKmers = new HashSet<>(bound);
        
        /* not connected, search from right */
        ArrayDeque<Kmer> rightPath = new ArrayDeque<>(bound);
        best = right;
        for (int depth=0; depth < bound; ++depth) {
            neighbors = best.getPredecessors(k, numHash, graph, minKmerCov);
//            best.predecessors = null; // clear cache
            
            if (neighbors.isEmpty()) {
                break;
            }
            else {
                if (neighbors.size() == 1) {
                    best = neighbors.peek();
                }
                else {
                    best = greedyExtendLeftOnce(graph, neighbors, lookahead);
                }
                
                if (best.equals(left)) {
                    return rightPath;
                }
                else {
                    if (rightPathKmers.contains(best)) {
                        return null;
                    }
                    else if (leftPathKmers.contains(best)) {
                        /* right path intersects the left path */
                        
                        if (isLowComplexity2(best.toString())) {
                            return null;
                        }
                        
                        rightPath.addFirst(best);

                        Iterator<Kmer> itr = leftPath.descendingIterator();
                        Kmer kmer;
                        while (itr.hasNext()) {
                            kmer = itr.next();
                            if (best.equals(kmer)) {
                                while (itr.hasNext()) {
                                    rightPath.addFirst(itr.next());
                                }

                                return rightPath;
                            }
                        }
                    }
                    else {
                        rightPathKmers.add(best);
                        rightPath.addFirst(best);
                    }
                }
            }
        }
        
        return null;
    }
    
    public static ArrayDeque<Kmer> getMaxCoveragePath(final BloomFilterDeBruijnGraph graph, 
                                                    final Kmer left, 
                                                    final Kmer right,
                                                    final int bound, 
                                                    final int lookahead, 
                                                    final BloomFilter bf) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        HashSet<Kmer> leftPathKmers = new HashSet<>(bound);
        
        /* extend right */
        ArrayDeque<Kmer> leftPath = new ArrayDeque<>(bound);
        Kmer best;
        ArrayDeque<Kmer> neighbors;
        
        best = left;

        for (int depth=0; depth < bound; ++depth) {
            neighbors = best.getSuccessors(k, numHash, graph, bf);
//            best.successors = null; // clear cache

            if (neighbors.isEmpty()) {
                break;
            }
            else {
                if (neighbors.size() == 1) {
                    best = neighbors.peek();
                }
                else {
                    best = greedyExtendRightOnce(graph, neighbors, lookahead, bf);
                }

                if (best.equals(right)) {
                    return leftPath;
                }
                else {
                    if (leftPathKmers.contains(best)) {
                        break;
                    }
                    else {
                        leftPathKmers.add(best);
                        leftPath.add(best);
                    }
                }
            }
        }
        
        HashSet<Kmer> rightPathKmers = new HashSet<>(bound);
        
        /* not connected, search from right */
        ArrayDeque<Kmer> rightPath = new ArrayDeque<>(bound);
        best = right;
        for (int depth=0; depth < bound; ++depth) {
            neighbors = best.getPredecessors(k, numHash, graph, bf);
//            best.predecessors = null; // clear cache
            
            if (neighbors.isEmpty()) {
                break;
            }
            else {
                if (neighbors.size() == 1) {
                    best = neighbors.peek();
                }
                else {
                    best = greedyExtendLeftOnce(graph, neighbors, lookahead, bf);
                }
                
                if (best.equals(left)) {
                    return rightPath;
                }
                else if (leftPathKmers.contains(best)) {
                    /* right path intersects the left path */
                    rightPath.addFirst(best);
                    
                    Iterator<Kmer> itr = leftPath.descendingIterator();
                    Kmer kmer;
                    while (itr.hasNext()) {
                        kmer = itr.next();
                        if (best.equals(kmer)) {
                            while (itr.hasNext()) {
                                rightPath.addFirst(itr.next());
                            }
                            
                            return rightPath;
                        }
                    }
                }
                else if (!rightPathKmers.contains(best)) {
                    rightPathKmers.add(best);
                    rightPath.addFirst(best);
                }
                else {
                    return null;
                }
            }
        }
        
        return null;
    }
        
    private static float[] getMinMedMax(float[] a) {
        int len = a.length;
        Arrays.sort(a);
        int halfLen = len/2;
        if (len % 2 == 0) {
            return new float[]{a[0], (a[halfLen-1] + a[halfLen])/2.0f, a[len-1]};
        }
        
        return new float[]{a[0], a[halfLen], a[len-1]};
    }
    
//    public static float getMedian(float[] arr) {
//        int len = arr.length;
//        float[] a = Arrays.copyOf(arr, len);
//        Arrays.sort(a);
//        int halfLen = len/2;
//        if (len % 2 == 0) {
//            return (a[halfLen-1] + a[halfLen])/2.0f;
//        }
//        
//        return a[halfLen];
//    }
    
    private static float getMedian(float[] arr) {
        int len = arr.length;
        Arrays.sort(arr);
        int halfLen = len/2;
        if (len % 2 == 0) {
            return (arr[halfLen-1] + arr[halfLen])/2.0f;
        }
        
        return arr[halfLen];
    }
    
    public static float getMinium(float[] arr) {
        float min = Float.MAX_VALUE;
        for (float c : arr) {
            if (c < min) {
                min = c;
            }
        }
        return min;
    }
        
    public static float getMedian(ArrayDeque<Float> arr) {
        int len = arr.size();
        ArrayList<Float> a = new ArrayList<>(arr);
        Collections.sort(a);
        int halfLen = len/2;
        if (len % 2 == 0) {
            return (a.get(halfLen-1) + a.get(halfLen))/2.0f;
        }
        
        return a.get(halfLen);
    }

    public static class LengthStats {
        public int min;
        public int q1;
        public int median;
        public int q3;
        public int max;
    }
    
    public static LengthStats getLengthStats(final int[] lengths) {
        int len = lengths.length;
        
        Arrays.sort(lengths);
        
        LengthStats stats = new LengthStats();
        int halfLen = len/2;
        int q1Index = len/4;
        int q3Index = halfLen+q1Index;
                
        stats.min = lengths[0];
        stats.max = lengths[len-1];
        
        if (len % 2 == 0) {
            stats.median = (lengths[halfLen-1] + lengths[halfLen])/2;
        }
        else {
            stats.median = lengths[halfLen];
        }
        
        if (len % 4 == 0) {
            stats.q1 = (lengths[q1Index-1] + lengths[q1Index])/2;
            stats.q3 = (lengths[q3Index-1] + lengths[q3Index])/2;
        }
        else {
            stats.q1 = lengths[q1Index];
            stats.q3 = lengths[q3Index];
        }
        
        return stats;
    }
    
    public static class CoverageStats {
        public float min;
        public float q1;
        public float median;
        public float q3;
        public float max;
        public float dropoff;
    }
    
    public static CoverageStats getCoverageStats(final Collection<Kmer> kmers, float maxCovGradient, int lookahead, boolean calcDropOff) {
        int len = kmers.size();
        float[] covs = new float[len];
        int i=0;
        for (Kmer kmer : kmers) {
            covs[i] = kmer.count;
            ++i;
        }
        
        Arrays.sort(covs);
                
        CoverageStats stats = new CoverageStats();
        int halfLen = len/2;
        int q1Index = len/4;
        int q3Index = halfLen+q1Index;
                
        stats.min = covs[0];
        stats.max = covs[len-1];
        
        if (len % 2 == 0) {
            stats.median = (covs[halfLen-1] + covs[halfLen])/2.0f;
        }
        else {
            stats.median = covs[halfLen];
        }
        
        if (len % 4 == 0) {
            stats.q1 = (covs[q1Index-1] + covs[q1Index])/2.0f;
            stats.q3 = (covs[q3Index-1] + covs[q3Index])/2.0f;
        }
        else {
            stats.q1 = covs[q1Index];
            stats.q3 = covs[q3Index];
        }
        
        stats.dropoff = 0;
        
        if (calcDropOff && len >= lookahead) {
            float last = covs[len-lookahead];
            for (i=len-lookahead-1; i>=0; --i) {
                float c = covs[i];
                if (c < last * maxCovGradient) {
                    stats.dropoff = last;
                    break;
                }
                last = c;
            }
        }
        
        return stats;
    }
    
    private static float rightGuidedMedianCoverage(BloomFilterDeBruijnGraph graph, String source, String guide, int guideLen) {
        if (guideLen == 0) {
            return graph.getCount(source);
        }
        
        float[] covs = new float[guideLen+1];
        covs[guideLen] = graph.getCount(source);
        
        String postfix = source.substring(1);
        String kmer;
        float count;
        for (int i=0; i<guideLen; ++i) {
            kmer = postfix + guide.charAt(i);
            count = graph.getCount(kmer);
            if (count > 0) {
                covs[i] = count;
                postfix = kmer.substring(1);
            }
            else {
                // not a valid sequence
                return 0;
            }
        }
        
        return getMedian(covs);
    }
    
    private static float leftGuidedMedianCoverage(BloomFilterDeBruijnGraph graph, String source, String guide, int guideLen) {
        if (guideLen == 0) {
            return graph.getCount(source);
        }
        
        float[] covs = new float[guideLen+1];
        covs[guideLen] = graph.getCount(source);
        
        int kMinus1 = graph.getK()-1;
        String prefix = source.substring(0,kMinus1);
        String kmer;
        float count;
        for (int i=guideLen-1; i>0; --i) {
            kmer = guide.charAt(i) + prefix;
            count = graph.getCount(kmer);
            if (count > 0) {
                covs[i] = count;
                prefix = kmer.substring(0,kMinus1);
            }
            else {
                // not a valid sequence
                return 0;
            }
        }
        
        return getMedian(covs);
    }
    
    public static ArrayDeque<Kmer> greedyExtendLeft(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead, int bound) {
        ArrayDeque<Kmer> extension = new ArrayDeque<>(bound);
        
        Kmer nextKmer = source;
        for (int i=0; i<bound; ++i) {
            nextKmer = greedyExtendLeftOnce(graph, nextKmer, lookahead);
            
            if (nextKmer == null) {
                break;
            }
            
            extension.addFirst(nextKmer);
        }
        
        return extension;
    }
    
    public static ArrayDeque<Kmer> greedyExtendLeftReversed(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead, int bound) {
        ArrayDeque<Kmer> extension = new ArrayDeque<>(bound);
        
        Kmer nextKmer = source;
        for (int i=0; i<bound; ++i) {
            nextKmer = greedyExtendLeftOnce(graph, nextKmer, lookahead);
            
            if (nextKmer == null) {
                break;
            }
            
            extension.addLast(nextKmer);
        }
        
        return extension;
    }
    
    public static ArrayDeque<Kmer> greedyExtendLeft(final BloomFilterDeBruijnGraph graph,
                                                    final Kmer source,
                                                    final int lookahead,
                                                    final int bound,
                                                    final BloomFilter bf) {
        ArrayDeque<Kmer> extension = new ArrayDeque<>(bound);
        
        Kmer nextKmer = source;
        for (int i=0; i<bound; ++i) {
            nextKmer = greedyExtendLeftOnce(graph, nextKmer, lookahead, bf);
            
            if (nextKmer == null) {
                break;
            }
            
            extension.addFirst(nextKmer);
        }
        
        return extension;
    }
    
    public static ArrayDeque<Kmer> greedyExtendRight(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead, int bound) {
        ArrayDeque<Kmer> extension = new ArrayDeque<>(bound);
        
        Kmer nextKmer = source;
        for (int i=0; i<bound; ++i) {
            nextKmer = greedyExtendRightOnce(graph, nextKmer, lookahead);
            
            if (nextKmer == null) {
                break;
            }
            
            extension.addLast(nextKmer);
        }
        
        return extension;
    }
    
    public static ArrayDeque<Kmer> greedyExtendRight(final BloomFilterDeBruijnGraph graph,
                                                    final Kmer source, 
                                                    final int lookahead, 
                                                    final int bound, 
                                                    final BloomFilter bf) {
        ArrayDeque<Kmer> extension = new ArrayDeque<>(bound);
        
        Kmer nextKmer = source;
        for (int i=0; i<bound; ++i) {
            nextKmer = greedyExtendRightOnce(graph, nextKmer, lookahead, bf);
            
            if (nextKmer == null) {
                break;
            }
            
            extension.addLast(nextKmer);
        }
        
        return extension;
    }
//    
//    public static String correctMismatches(String seq, BloomFilterDeBruijnGraph graph, int lookahead, int mismatchesAllowed) {
//        int numCorrected = 0;
//        int seqLen = seq.length();
//        int k = graph.getK();
//        
//        if (seqLen < k) {
//            // no correction
//            return seq;
//        }
//        
//        int numKmers = seqLen-k+1;
//        StringBuilder sb = new StringBuilder(seq);
//        
//        float bestCov, cov;
//        String kmer, guide;
//        ArrayDeque<String> variants;
//        char bestBase = 'N';
//        int guideEnd = -1;
//        int guideLen = -1;
//        
//        // correct from start
//        for (int i=0; i<numKmers; ++i) {
//            int j = i+k;
//            kmer = sb.substring(i, j);
//            variants = graph.getRightVariants(kmer);
//            if (!variants.isEmpty()) {
//                guideEnd = Math.min(j+lookahead, seqLen);
//                guide = sb.substring(j, guideEnd);
//                guideLen = guideEnd - j;
//                bestCov = 0;
//                
//                if (graph.contains(kmer)) {
//                    bestCov = rightGuidedMedianCoverage(graph, kmer, guide, guideLen);
//                }
//                
//                boolean corrected = false;
//                for (String v : variants) {
//                    cov = rightGuidedMedianCoverage(graph, v, guide, guideLen);
//                    if (cov > bestCov) {
//                        bestCov = cov;
//                        bestBase = v.charAt(k-1);
//                        corrected = true;
//                    }
//                }
//                
//                if (corrected) {
//                    if (++numCorrected > mismatchesAllowed) {
//                        // too many mismatches
//                        return seq;
//                    }
//                    else {
//                        sb.setCharAt(j-1, bestBase);
//                    }
//                    
//                    //i += lookahead;
//                }
//            }
//        }
//        
//        // correct from end
//        for (int i=seqLen-k; i>=0; --i) {
//            kmer = sb.substring(i, i+k);
//            variants = graph.getLeftVariants(kmer);
//            if (!variants.isEmpty()) {
//                guideEnd = Math.max(0, i-lookahead);
//                guideLen = i - guideEnd;
//                guide = sb.substring(guideEnd, i);
//                bestCov = 0;
//                
//                if (graph.contains(kmer)) {
//                    bestCov = leftGuidedMedianCoverage(graph, kmer, guide, guideLen);
//                }
//                
//                boolean corrected = false;
//                for (String v : variants) {
//                    cov = leftGuidedMedianCoverage(graph, v, guide, guideLen);
//                    if (cov > bestCov) {
//                        bestCov = cov;
//                        bestBase = v.charAt(0);
//                        corrected = true;
//                    }
//                }
//                
//                if (corrected) {                    
//                    if (++numCorrected > mismatchesAllowed) {
//                        // too many mismatches
//                        return seq;
//                    }
//                    else {
//                        sb.setCharAt(i, bestBase);
//                    }
//                    
//                    //i -= lookahead;
//                }
//            }
//        }
//        
//        if (numCorrected == 0) {
//            return seq;
//        }
//        
//        String seq2 = sb.toString();
//        if (graph.isValidSeq(seq2)) {
//            return seq2;
//        }
//        
//        return seq;
//    }
    
    private static int getHighCoverageWindowIndex(ArrayList<Kmer> kmers, int start, int window, int halfWindow, float covThreshold) {
        int end = kmers.size() - window;
        
        for (int i=start; i<end; ++i) {
            if (kmers.get(i).count >= covThreshold) {
                boolean found = true;
                for (int j=i+1; j<window; ++j) {
                    if (kmers.get(j).count < covThreshold) {
                        found = false;
                        i = j;
                    }
                }
                
                if (found) {
                    return i + halfWindow;
                }
            }
        }
        
        return start;
    }
    
    private static ArrayList<Kmer> correctErrorHelper3(final ArrayList<Kmer> kmers,
                                                    final BloomFilterDeBruijnGraph graph, 
                                                    final int lookahead,
                                                    final int maxIndelSize,
                                                    final float maxCovGradient,
                                                    final float percentIdentity) {
        final int numKmers = kmers.size();
        final int windowSize = 3;
        
        if (numKmers > windowSize) {
            final int k = graph.getK();
            final int numHash = graph.getMaxNumHash();
            final int snvBubbleLength = k;
            
            ArrayList<Kmer> kmers2 = new ArrayList<>(numKmers + maxIndelSize);

            float[] windowCounts = new float[windowSize];
            for (int i=0; i<windowSize; ++i) {
                windowCounts[i] = kmers.get(i).count;
            }
            
            FloatRollingWindow window = new FloatRollingWindow(windowCounts);
            float lastMin = window.getMin();
            float currentMin;
            int lastGoodKmerIndex = -1;
            
            for (int i=windowSize; i<numKmers; ++i) {
                
                window.roll(kmers.get(i).count);
                currentMin = window.getMin();
                
                if (currentMin <= lastMin) {
                    if (currentMin <= lastMin * maxCovGradient) {
                        // ... ---___

                        for (int j=lastGoodKmerIndex+1; j<i; ++j) {
                            kmers2.add(kmers.get(j));
                        }
                        
                        lastGoodKmerIndex = i;
                    }
                }
                else {
                    if (lastMin <= currentMin * maxCovGradient) {
                        // ... ___---
                        
                        if (lastGoodKmerIndex < 0) {
                            // ______---
                            
                            // correct left edge
                            // TODO
                        }
                        else {
                            // ---______---
                            
                            // bridge gap
                            // TODO
                        }
                    }
                }
            }
        }
        
        return null;
    }
    
    private static ArrayList<Kmer> correctErrorHelper2(final ArrayList<Kmer> kmers,
                                                    final BloomFilterDeBruijnGraph graph, 
                                                    final int lookahead,
                                                    final int maxIndelSize,
                                                    final float covThreshold,
                                                    final float percentIdentity) {
        
        /**@TODO Need to be debugged*/
        
        final int numKmers = kmers.size();
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int snvBubbleLength = k;
        final int minIslandSize = 3;
        
        ArrayList<Kmer> kmers2 = new ArrayList<>(numKmers + maxIndelSize);
        
        boolean corrected = false;
                
        int start = -1; // left edge of kmer island
        int end = -1;   // right edge of kmer island
        
        for (int i=0; i<numKmers; ++i) {
            if (start < 0) {
                // the left edge of kmer island has not been defined yet
                boolean enoughKmers = true;
                for (int j=i; j<i+minIslandSize; ++j) {
                    if (kmers.get(j).count < covThreshold) {
                        enoughKmers = false;
                        break;
                    }
                }
                
                if (enoughKmers) {
                    start = i;
                    
                    // attempt to correct bad kmers to the left of island
                    
                    if (end < 0) {
                        if (start > 0) {
                            // correct left edge

                            ArrayDeque<Kmer> leftVars = kmers.get(start-1).getLeftVariants(k, numHash, graph);
                            
                            if (leftVars.isEmpty()) {
                                // no branches found
                                for (int j=0; j<start; ++j) {
                                    kmers2.add(kmers.get(j));
                                }
                            }
                            else if (start >= lookahead) {
                                // calculate median cov of edge kmers
                                float[] tipCovs = new float[start];
                                for (int j=0; j<start; ++j) {
                                    tipCovs[j] = kmers.get(j).count;
                                }
                                float tipMedCov = getMedian(tipCovs);

                                ArrayDeque<Kmer> greedyTipKmers = greedyExtendLeft(graph, kmers.get(start), lookahead, start);
                                if (greedyTipKmers.size() == start && getMedianKmerCoverage(greedyTipKmers) > tipMedCov) {
                                    if (getPercentIdentity(graph.assemble(greedyTipKmers), graph.assemble(kmers, 0, start)) >= percentIdentity){
                                        corrected = true;
                                        kmers2.addAll(greedyTipKmers);
                                    }
                                    else if (!kmers.get(0).hasPredecessors(k, numHash, graph) && start < k) {
                                        // this is blunt end in graph
                                        // do not add its kmers
                                        corrected = true;
                                    }
                                    else {
                                        // use original sequence
                                        for (int j=0; j<start; ++j) {
                                            kmers2.add(kmers.get(j));
                                        }
                                    }
                                }
                                else {
                                    // use original sequence
                                    for (int j=0; j<start; ++j) {
                                        kmers2.add(kmers.get(j));
                                    }
                                }
                            }
                            else {
                                corrected = true;
                            }                        
                        }
                    }
                    else {
                        // correct bubble
                        int bubbleLength = start - end - 1;
                        
                        if (bubbleLength == snvBubbleLength) {
                            // a SNV bubble
                            ArrayList<Kmer> bestKmers = null;
                            float bestCov = Float.MIN_VALUE;
                            
                            String prefix = graph.getSuffix(kmers.get(end).toString());
                            String suffix = graph.getPrefix(kmers.get(start).toString());
                            
                            for (char n : NUCLEOTIDES) {
                                ArrayList<Kmer> testKmers = graph.getKmers(prefix + n + suffix);
                                float[] m3 = getMinMedMaxKmerCoverage(testKmers);

                                float medCov = m3[1];
                                if (m3[0] > 0 && medCov > bestCov) {
                                    bestCov = medCov;
                                    bestKmers = testKmers;
                                }
                            }

                            if (bestKmers != null & bestCov > 0) {
                                // fill with best kmers
                                kmers2.addAll(bestKmers);
                                corrected = true;
                            }
                            else {
                                // fill with original kmers
                                for (int j=end+1; j<start; ++j) {
                                    kmers2.add(kmers.get(j));
                                }                            
                            }
                        }
                        else {
                            // non-SNV bubble
                            
                            ArrayDeque<Kmer> path = getMaxCoveragePath(graph, kmers2.get(end+1), kmers.get(start-1), bubbleLength + maxIndelSize, lookahead, 1);
                            if (path == null) {
                                // fill with original sequence
                                for (int j=end+1; j<start; ++j) {
                                    kmers2.add(kmers.get(j));
                                }
                            }
                            else {
                                int altPathLen = path.size();

                                if (bubbleLength-maxIndelSize <= altPathLen && 
                                        altPathLen <= bubbleLength+maxIndelSize && 
                                        (altPathLen <= k+maxIndelSize ||
                                            getPercentIdentity(graph.assemble(path), graph.assemble(kmers, end+1, start)) >= percentIdentity)) {

                                        kmers2.addAll(path);
                                        corrected = true;
                                }
                                else {
                                    // fill with original sequence
                                    for (int j=end+1; j<start; ++j) {
                                        kmers2.add(kmers.get(j));
                                    }
                                }
                            }
                        }
                        
                        end = i;
                    }
                    
                    kmers2.add(kmers.get(i));
                }
            }
            else if (kmers.get(i).count >= covThreshold) {
                end = i;
                kmers2.add(kmers.get(i));
            }
            else {
                start = -1;
            }
        }
        
        if (end < numKmers-1) {
            // correct right edge
            
            int rightTipLen = numKmers-1-end;
            
            int firstBadKmerIndex = end + 1;
            
            ArrayDeque<Kmer> rightVars = kmers.get(firstBadKmerIndex).getRightVariants(k, numHash, graph);
            if (rightVars.isEmpty()) {
                for (int j=firstBadKmerIndex; j<numKmers; ++j) {
                    kmers2.add(kmers.get(j));
                }
            }
            else if (rightTipLen >= lookahead) {
                // calculate median cov of edge kmers
                float[] tipCovs = new float[rightTipLen];
                for (int j=0; j<rightTipLen; ++j) {
                    tipCovs[j] = kmers.get(firstBadKmerIndex+j).count;
                }
                float tipMedCov = getMedian(tipCovs);
                
                ArrayDeque<Kmer> greedyTipKmers = greedyExtendRight(graph, kmers.get(firstBadKmerIndex), lookahead, rightTipLen);
                if (greedyTipKmers.size() == rightTipLen && getMedianKmerCoverage(greedyTipKmers) > tipMedCov) {
                    if (getPercentIdentity(graph.assemble(greedyTipKmers), graph.assemble(kmers, firstBadKmerIndex, numKmers)) >= percentIdentity){
                        corrected = true;
                        kmers2.addAll(greedyTipKmers);
                    }
                    else if (!kmers.get(numKmers-1).hasSuccessors(k, numHash, graph) && rightTipLen < k) {
                        // this is blunt end in graph
                        // do not add its kmers
                        corrected = true;
                    }
                    else {
                        // use original sequence
                        for (int j=firstBadKmerIndex; j<numKmers; ++j) {
                            kmers2.add(kmers.get(j));
                        }
                    }
                }
                else {
                    // use original sequence
                    for (int j=firstBadKmerIndex; j<numKmers; ++j) {
                        kmers2.add(kmers.get(j));
                    }
                }
            }
            else {
                corrected = true;
            }
            
        }
        
        corrected = correctMismatches(kmers2, graph, covThreshold, 1) || corrected;
        
        if (corrected) {
            return kmers2;
        }
        
        return null;
    }
    
    private static void updateSketch(long[] sketch, int sketchSize, long newVal) {
        int i = sketchSize - 1;
        if (newVal >= sketch[i]) {
            return;
        }
        
        // find index to insert new hash value
        for (--i; i>=0; --i) {
            long val = sketch[i];
            
            if (newVal == val) {
                return;
            }
            else if (newVal > val) {
                break;
            }
        }
        
        ++i;
        
        // shift larger hash values
        System.arraycopy(sketch, i, sketch, i+1, sketchSize-i-1);
        
        // insert the new hash value
        sketch[i] = newVal;
    }
    
    public static long[] getMinimizersWithCompressedHomoPolymers(String seq, int k, NTHashIterator itr, int windowSize) {
        String chpSeq = compressHomoPolymers(seq);
        return getMinimizers(chpSeq, getNumKmers(chpSeq, k), itr, windowSize);
    }

    public static TreeSet<Long> getMinimizersSetWithCompressedHomoPolymers(String seq, int k, NTHashIterator itr, int windowSize) {
        String chpSeq = compressHomoPolymers(seq);
        return getMinimizersSet(chpSeq, getNumKmers(chpSeq, k), itr, windowSize);
    }
    
    public static long[] getMinimizers(String seq, int numKmers, NTHashIterator itr, int windowSize) {        
        TreeSet<Long> minimizers = getMinimizersSet(seq, numKmers, itr, windowSize);
        
        long[] minimizersArr = new long[minimizers.size()];
        int i=0;
        for (Long m : minimizers) {
            minimizersArr[i++] = m;
        }
        
        return minimizersArr;
    }
    
    public static TreeSet<Long> getMinimizersSet(String seq, int numKmers, NTHashIterator itr, int windowSize) {
        itr.start(seq);
        long[] hvals = itr.hVals;

        if (numKmers <= windowSize) {
            long minimizer = hvals[0];
            while(itr.hasNext()) {
                itr.next();
                long h = hvals[0];
                if (h < minimizer) {
                    minimizer = h;
                }
            }
            
            TreeSet<Long> set = new TreeSet<>();
            set.add(minimizer);
            return set;
        }
        
        TreeSet<Long> minimizers = new TreeSet<>();
        
        // find minimizer in the first window
        ArrayDeque<Long> window = new ArrayDeque<>(windowSize);
        
        itr.next();
        long minimizer = hvals[0];
        int minmizerWindowPos = 0;
        window.add(minimizer);
        
        for (int i=1; i<windowSize; ++i) {
            itr.next();
            long h = hvals[0];
            window.add(h);
            
            if (h < minimizer) {
                minimizer = h;
                minmizerWindowPos = i;
            }
        }
        minimizers.add(minimizer);
        
        // find minimizers as window shifts
        while (itr.hasNext()) {
            window.removeFirst();
            
            itr.next();
            long h = hvals[0];
            window.add(h);
            
            if (--minmizerWindowPos < 0) {
                Iterator<Long> wItr = window.iterator();
                
                minimizer = wItr.next();
                minmizerWindowPos = 0;
                
                for(int i=1; i<windowSize; ++i) {
                    h = wItr.next();
                    if (h < minimizer) {
                        minimizer = h;
                        minmizerWindowPos = i;
                    }
                }
                
                minimizers.add(minimizer);
            }
            else if (h < minimizer) {
                minimizer = h;
                minmizerWindowPos = windowSize-1;
                minimizers.add(minimizer);
            }
        }
        
        return minimizers;
    }
    
    public static long[] getMinimizers(String seq, int numKmers, NTHashIterator itr, int windowSize, BloomFilterDeBruijnGraph graph, float minCoverage) {        
        TreeSet<Long> minimizers = getMinimizersSet(seq, numKmers, itr, windowSize, graph, minCoverage);
        
        long[] minimizersArr = new long[minimizers.size()];
        int i=0;
        for (Long m : minimizers) {
            minimizersArr[i++] = m;
        }
        
        return minimizersArr;
    }
    
    public static TreeSet<Long> getMinimizersSet(String seq, int numKmers, NTHashIterator itr, int windowSize, BloomFilterDeBruijnGraph graph, float minCoverage) {        
        itr.start(seq);
        long[] hvals = itr.hVals;

        if (numKmers <= windowSize) {
            long minimizer = hvals[0];
            boolean hasMinimizer = graph.getCount(hvals) >= minCoverage;
            
            while(itr.hasNext()) {
                itr.next();
                long h = hvals[0];
                
                if (!hasMinimizer) {
                    if (graph.getCount(hvals) >= minCoverage) {
                        minimizer = h;
                        hasMinimizer = true;
                    }
                }
                else if (h < minimizer) {
                    if (graph.getCount(hvals) >= minCoverage) {
                        minimizer = h;
                    }
                }
            }
            
            TreeSet<Long> set = new TreeSet<>();
            if (hasMinimizer) {
                set.add(minimizer);
            }
            return set;
        }
        
        TreeSet<Long> minimizers = new TreeSet<>();
        ArrayDeque<Long> window = new ArrayDeque<>(windowSize);
        ArrayDeque<Float> covWindow = new ArrayDeque<>(windowSize);
        
        // find minimizer in the first window
        itr.next();
        int minimizerWindowPos = 0;
        long minimizer = hvals[0];
        window.add(minimizer);
        float c = graph.getCount(hvals);
        covWindow.add(c);
        boolean hasMinimizer = c >= minCoverage;
        
        for (int i=1; i<windowSize; ++i) {
            itr.next();
            
            long h = hvals[0];
            window.add(h);
            
            c = graph.getCount(hvals);
            covWindow.add(c);
            
            if (c >= minCoverage) {
                if (!hasMinimizer) {
                    minimizer = h;
                    minimizerWindowPos = i;
                    hasMinimizer = true;
                }
                else if (h < minimizer) {
                    minimizer = h;
                    minimizerWindowPos = i;
                }
            }
        }
        
        if (hasMinimizer) {
            minimizers.add(minimizer);
        }
        
        // find minimizers as window shifts
        while (itr.hasNext()) {
            window.removeFirst();
            covWindow.removeFirst();
            --minimizerWindowPos;
            
            itr.next();
            
            long h = hvals[0];
            window.add(h);
            
            c = graph.getCount(hvals);
            covWindow.add(c);
            
            if (c >= minCoverage) {
                if (!hasMinimizer) {
                    minimizer = h;
                    minimizerWindowPos = windowSize-1;
                    minimizers.add(minimizer);
                    hasMinimizer = true;
                }
                else if (minimizerWindowPos < 0) {
                    Iterator<Long> wItr = window.iterator();
                    Iterator<Float> cItr = covWindow.iterator();
                    
                    minimizerWindowPos = 0;
                    minimizer = wItr.next();
                    c = cItr.next();
                    hasMinimizer = c >= minCoverage;
                    
                    for(int i=1; i<windowSize; ++i) {
                        h = wItr.next();
                        c = cItr.next();
                        
                        if (c >= minCoverage) {
                            if (!hasMinimizer) {
                                minimizer = h;
                                minimizerWindowPos = i;
                                hasMinimizer = true;
                            }
                            else if (h < minimizer) {
                                minimizer = h;
                                minimizerWindowPos = i;
                            }
                        }
                    }
                    
                    if (hasMinimizer) {
                        minimizers.add(minimizer);
                    }
                }
                else if (h < minimizer) {
                    minimizer = h;
                    minimizerWindowPos = windowSize-1;
                    minimizers.add(minimizer);
                }
            }
        }
        
        return minimizers;
    }
    
    public static long[] getAscendingHashValuesWithCompressedHomoPolymers(String seq, NTHashIterator itr, int k) {
        TreeSet<Long> hashValSet = getHashValuesSetWithCompressedHomoPolymers(seq, itr, k);
        
        int numVals = hashValSet.size();
        long[] result = new long[numVals];
        int i=0;
        for (Long h : hashValSet) {
            result[i++] = h;
        }
        
        return result;
    }
 
    public static TreeSet<Long> getHashValuesSetWithCompressedHomoPolymers(String seq, NTHashIterator itr, int k) {
        seq = compressHomoPolymers(seq);
        int numKmers = getNumKmers(seq, k);
        
        TreeSet<Long> hashValSet = new TreeSet<>();
        
        itr.start(seq);
        long[] hVals = itr.hVals;
        for (int i=0; i<numKmers; ++i) {
            itr.next();
            hashValSet.add(hVals[0]);
        }
        
        return hashValSet;
    }
    
    public static long[] getAscendingHashValues(String seq, NTHashIterator itr, BloomFilterDeBruijnGraph graph, int numKmers, float minCoverage) {
        TreeSet<Long> hashValSet = getHashValuesSet(seq, itr, graph, numKmers, minCoverage);
        
        int numVals = hashValSet.size();
        long[] result = new long[numVals];
        int i=0;
        for (Long h : hashValSet) {
            result[i++] = h;
        }
        
        return result;
    }
    
    public static TreeSet<Long> getHashValuesSet(String seq, NTHashIterator itr, BloomFilterDeBruijnGraph graph, int numKmers, float minCoverage) {
        TreeSet<Long> hashValSet = new TreeSet<>();
        
        itr.start(seq);
        long[] hVals = itr.hVals;
        for (int i=0; i<numKmers; ++i) {
            itr.next();
            if (graph.getCount(hVals) >= minCoverage) {
                hashValSet.add(hVals[0]);
            }
        }
        
        return hashValSet;
    }
    
    public static long[] combineSketches(long[]... sketches) {
        TreeSet<Long> hashValSet = new TreeSet<>();
        for (long[] sketch : sketches) {
            for (long val : sketch) {
                hashValSet.add(val);
            }
        }
        
        int numVals = hashValSet.size();
        long[] result = new long[numVals];
        int i=0;
        for (Long h : hashValSet) {
            result[i++] = h;
        }
        
        return result;
    }
    
    public static long[] combineSketches(Collection<long[]> sketches) {
        TreeSet<Long> hashValSet = new TreeSet<>();
        for (long[] sketch : sketches) {
            for (long val : sketch) {
                hashValSet.add(val);
            }
        }
        
        int numVals = hashValSet.size();
        long[] result = new long[numVals];
        int i=0;
        for (Long h : hashValSet) {
            result[i++] = h;
        }
        
        return result;
    }
    
    public static long[] getBottomSketch(long[] sortedHashVals, int sketchSize) {
        long[] sketch = new long[sketchSize];
        System.arraycopy(sortedHashVals, 0, sketch, 0, sketchSize);
        return sketch;
    }
        
    public static int getNumIntersection(long[] sketch1, long[] sketch2) {
        int sketchSize1 = sketch1.length;
        int sketchSize2 = sketch2.length;
        int intersectionSize = 0;
        int j = 0;
        
        long hash1, hash2;
        for (int i=0; i<sketchSize1 && j<sketchSize2; ++i) {
            hash1 = sketch1[i];
            hash2 = sketch2[j];
            
            if (hash1 == hash2) {
                ++intersectionSize;
                ++j;
            }
            else if (hash1 > hash2) {
                for (++j; j<sketchSize2; ++j) {
                    hash2 = sketch2[j];
                    if (hash1 == hash2) {
                        ++intersectionSize;
                        ++j;
                        break;
                    }
                    else if (hash1 < hash2) {
                        break;
                    }
                }
            }
        }
        
        return intersectionSize;
    }
    
    public static int getNumIntersection(TreeSet<Long> sketch1, TreeSet<Long> sketch2) {
        Iterator<Long> itr1 = sketch1.iterator();
        Iterator<Long> itr2 = sketch2.iterator();
        int intersectionSize = 0;
        
        if (!itr1.hasNext() || !itr2.hasNext()) {
            return 0;
        }
        
        long hash1 = itr1.next();
        long hash2 = itr2.next();
        
        while (true) {
            if (hash1 == hash2) {
                ++intersectionSize;
                
                if (itr1.hasNext()) {
                    hash1 = itr1.next();
                }
                else {
                    break;
                }
                
                if (itr2.hasNext()) {
                    hash2 = itr2.next();
                }
                else {
                    break;
                }
            }
            else if (hash1 > hash2) {
                if (itr2.hasNext()) {
                    hash2 = itr2.next();
                }
                else {
                    break;
                }
            }
            else if (hash1 < hash2) {
                if (itr1.hasNext()) {
                    hash1 = itr1.next();
                }
                else {
                    break;
                }
            }
        }
        
        return intersectionSize;
    }
    
    /**
     * Method to return the jaccard similarity coefficient of 2 bottom sketches.
     * The 2 sketches must have the same size.
     * @param sketch1 sorted bottom sketch from set 1
     * @param sketch2 sorted bottom sketch from set 2
     * @return jaccard similarity coefficient
     */
    public static float getResemblance(long[] sketch1, long[] sketch2) {
        int intersectionSize = getNumIntersection(sketch1, sketch2);
        
        return getNumIntersection(sketch1, sketch2) / (float) (sketch1.length + sketch2.length - intersectionSize);
    }
    
    public static float getContainment(long[] test, long[] target) {        
        return getNumIntersection(test, target) / (float) (target.length);
    }
    
    public static float getAdjustedMeanKmerCoverage(ArrayList<Kmer> kmers, float maxCovGradient, int lookahead) {
        CoverageStats covStat = getCoverageStats(kmers, maxCovGradient, lookahead, false);
        float maxCov = covStat.median + (covStat.q3 - covStat.q1) * 1.5f;
        
        int numKmers = 0;
        float total = 0;
        
        for (Kmer kmer : kmers) {
            float c = kmer.count;
            if (c <= maxCov) {
                total += c;
                ++numKmers;
            }
        }
        
        return total/numKmers; 
    }
    
    public static float getMinMedianKmerCoverage(ArrayList<Kmer> kmers, int windowSize) {
        int numKmers = kmers.size();
        
        if (numKmers <= windowSize) {
            return getMedianKmerCoverage(kmers);
        }
        
        float minCov = Float.MAX_VALUE;
        
        float[] window = new float[windowSize];
        for (int i=(numKmers % windowSize)/2; i+windowSize<=numKmers; i+=windowSize) {
            for (int j=0; j<windowSize; ++j) {
                window[j] = kmers.get(i+j).count;
            }
            
            float m = getMedian(window);
            if (m < minCov) {
                minCov = m;
            }
        }
                
        return minCov;
    }
    
    public static ArrayList<Kmer> correctLongSequence(ArrayList<Kmer> kmers, 
                                                    BloomFilterDeBruijnGraph graph, 
                                                    int maxErrCorrItr, 
                                                    float maxCovGradient, 
                                                    int lookahead, 
                                                    int maxIndelSize, 
                                                    float percentIdentity, 
                                                    int minKmerCov,
                                                    int minNumSolidKmers,
                                                    boolean trimLowCovEdges) {

        int numNeeded = minNumSolidKmers;
        for (Kmer kmer : kmers) {
            if (kmer.count >= minKmerCov) {
                if (--numNeeded <= 0) {
                    break;
                }
            }
        }

        if (numNeeded <= 0) {
            // use each window's median kmer coverage as threshold
            int windowSize = 500;
            int shift = windowSize/2;
            
            // trim head and tail
            ArrayList<Kmer> testKmers = trimLowCovEdges ?
                                        trimLowCoverageEdgeKmers(kmers, graph, lookahead, minKmerCov) :
                                        null;
            
            if (testKmers == null) {
                testKmers = kmers;
            }
            
            boolean toShift = false;
            for (int itr=0; itr<maxErrCorrItr; ++itr) {
                boolean corrected = false;
                
                int numKmers = testKmers.size();
                ArrayList<Kmer> correctedKmers = new ArrayList<>();

                int end = 0;
                for (int i=0; i<numKmers; i=end) {
                    if (i == 0 && toShift) { 
                        end = Math.min(i + shift + windowSize, numKmers);
                    }
                    else {
                        end = Math.min(i + windowSize, numKmers);
                    }
                    
                    if (end + shift >= numKmers) {
                        end = numKmers;
                    }

                    ArrayList<Kmer> window = new ArrayList<>(testKmers.subList(i, end));
                    ArrayList<Kmer> windowCorrected = null;

                    if (end-i >= lookahead) {
                        CoverageStats covStat = getCoverageStats(window, maxCovGradient, lookahead, true);
                        float threshold = covStat.median;
                        if (covStat.dropoff > 0) {
                            threshold = Math.max(covStat.dropoff, minKmerCov);
                        }

                        windowCorrected = correctInternalErrors(window,
                                                                    graph, 
                                                                    lookahead,
                                                                    maxIndelSize,
                                                                    threshold,
                                                                    percentIdentity,
                                                                    minKmerCov);
                    }
                    
                    if (windowCorrected == null) {
                        correctedKmers.addAll(window);
                    }
                    else {
                        corrected = true;
                        correctedKmers.addAll(windowCorrected);
                    }
                }
                
                if (!corrected) {
                    break;
                }

                testKmers = correctedKmers;
                toShift = !toShift;
            }
            
            int numKmers = testKmers.size();
            if (trimLowCovEdges && numKmers > lookahead) {
                ArrayList<Kmer> window = new ArrayList<>(testKmers.subList(0, Math.min(windowSize, numKmers)));
                CoverageStats headCovStat = getCoverageStats(window, maxCovGradient, lookahead, true);

                window = new ArrayList<>(testKmers.subList(Math.max(0, numKmers-windowSize), numKmers));
                CoverageStats tailCovStat = getCoverageStats(window, maxCovGradient, lookahead, true);

                if (headCovStat.dropoff > 0 || tailCovStat.dropoff > 0) {
                    ArrayList<Kmer> correctedKmers = correctEdgeErrors(testKmers,
                                                                        graph,
                                                                        lookahead,
                                                                        Math.max(minKmerCov, headCovStat.dropoff),
                                                                        Math.max(minKmerCov, tailCovStat.dropoff),
                                                                        percentIdentity,
                                                                        minKmerCov);
                    if (correctedKmers != null) {
                        testKmers = correctedKmers;
                    }
                }
            }
            
            return testKmers;
        }
        
        return null;
    }
    
    public static ArrayList<Kmer> trimLowCoverageEdgeKmers(ArrayList<Kmer> kmers,
                                                        BloomFilterDeBruijnGraph graph, 
                                                        int lookahead,
                                                        float threshold) {
        int k = graph.getK();
        int numKmers = kmers.size();
        int headIndex = 0;
        int tailIndex = numKmers;
        int window = k*5;
        float windowThreshold = (float)5/(float)k;

        for (int i=0; i<numKmers; ++i) {
            Kmer kmer = kmers.get(i);
            if (kmer.count >= threshold) {
                int end = Math.min(i+lookahead, numKmers);
                if (areKmerCoverageAboveThreshold(kmers, i+1, end, threshold)) {
                    if (end+window > numKmers || (float)getNumKmersAboveCoverageThreshold(kmers, end, end+window, threshold)/(float)window > windowThreshold) {
                        headIndex = i;
                        break;
                    }
                }
            }
        }

        for (int i=numKmers-1; i>headIndex; --i) {
            Kmer kmer = kmers.get(i);
            if (kmer.count >= threshold) {
                int start = Math.max(0, i-lookahead);
                if (areKmerCoverageAboveThreshold(kmers, start, i, threshold)) {
                    if (start-window < 0 || (float)getNumKmersAboveCoverageThreshold(kmers, start-window, start, threshold)/(float)window > windowThreshold) {
                        tailIndex = i+1;
                        break;
                    }
                }
            }
        }
        
        if (headIndex > 0 || tailIndex < numKmers) {
            if (tailIndex - headIndex > 1) {
                // trimming was done
                return new ArrayList<>(kmers.subList(headIndex, tailIndex));
            }
            else {
                // sequence was trimmed to nothing
                return new ArrayList<>();
            }
        }
        
        // no trimming was done
        return null;
    }
    
    public static ArrayList<Kmer> correctEdgeErrors(ArrayList<Kmer> kmers,
                                                    BloomFilterDeBruijnGraph graph,
                                                    int lookahead,
                                                    float headCovThreshold,
                                                    float tailCovThreshold,
                                                    float percentIdentity,
                                                    float minKmerCov) {
        int numKmers = kmers.size();
        int headIndex = 0;
        int tailIndex = numKmers-1;

        if (headCovThreshold > 0) {
            for (int i=0; i<numKmers; ++i) {
                Kmer kmer = kmers.get(i);
                if (kmer.count >= headCovThreshold) {
                    headIndex = i;
                    break;
                }
            }
        }
        
        if (tailCovThreshold > 0) {
            for (int i=numKmers-1; i>=0; --i) {
                Kmer kmer = kmers.get(i);
                if (kmer.count >= tailCovThreshold) {
                    tailIndex = i;
                    break;
                }
            }
        }
        
        ArrayDeque<Kmer> newHead = null;
        boolean trimHead = false;
        if (headIndex > 0) {
            int k = graph.getK();
            int numHash = graph.getMaxNumHash();
            
            Kmer original = kmers.get(headIndex-1);
            ArrayList<Kmer> path = new ArrayList<>(kmers.subList(0, headIndex));
            int pathLen = path.size();
            float pathCov = getMedianKmerCoverage(path);
            
            float altPathCov = 0;
            for (Kmer var : original.getLeftVariants(k, numHash, graph, minKmerCov)) {
                ArrayDeque<Kmer> varPath = greedyExtendLeft(graph, var, lookahead, pathLen-1);
                varPath.addLast(var);
                if (varPath.size() == pathLen) {
                    float cov = getMedianKmerCoverage(varPath);
                    if (cov > pathCov && cov > altPathCov && getPercentIdentity(graph.assemble(varPath), graph.assemble(kmers, 0, headIndex)) >= percentIdentity) {
                        newHead = varPath;
                    }
                }
            }
            
            if (newHead == null && headIndex < k) {
                trimHead = true;
            }
        }
        
        ArrayDeque<Kmer> newTail = null;
        boolean trimTail = false;
        if (tailIndex < numKmers-1) {
            int k = graph.getK();
            int numHash = graph.getMaxNumHash();
            
            Kmer original = kmers.get(tailIndex+1);
            ArrayList<Kmer> path = new ArrayList<>(kmers.subList(tailIndex+1, numKmers));
            int pathLen = path.size();
            float pathCov = getMedianKmerCoverage(path);
            
            float altPathCov = 0;
            for (Kmer var : original.getRightVariants(k, numHash, graph, minKmerCov)) {
                ArrayDeque<Kmer> varPath = greedyExtendRight(graph, var, lookahead, pathLen-1);
                varPath.addFirst(var);
                if (varPath.size() == pathLen) {
                    float cov = getMedianKmerCoverage(varPath);
                    if (cov > pathCov && cov > altPathCov &&
                            getPercentIdentity(graph.assemble(varPath), graph.assemble(kmers, tailIndex+1, numKmers)) >= percentIdentity) {
                        newTail = varPath;
                    }
                }
            }
            
            if (newTail == null && tailIndex >= numKmers - k) {
                trimTail = true;
            }
        }
        
        if (newHead != null || newTail != null || trimHead || trimTail) {
            ArrayList<Kmer> newKmers = new ArrayList<>(numKmers);
            
            if (trimHead) {
                newKmers.addAll(kmers.subList(headIndex, Math.min(tailIndex+1, numKmers)));
            }
            else if (newHead == null) {
                newKmers.addAll(kmers.subList(0, tailIndex+1));
            }
            else {
                newKmers.addAll(newHead);
                newKmers.addAll(kmers.subList(headIndex, Math.min(tailIndex+1, numKmers)));
            }
            
            if (!trimTail) {
                if (newTail == null) {
                    if (tailIndex+1 < numKmers) {
                        newKmers.addAll(kmers.subList(tailIndex+1, numKmers));
                    }
                }
                else {
                    newKmers.addAll(newTail);
                }
            }
            
            return newKmers;
        }
        
        return null;
    }
    
    public static ArrayList<Kmer> correctInternalErrors(ArrayList<Kmer> kmers,
                                                    BloomFilterDeBruijnGraph graph, 
                                                    int lookahead,
                                                    int maxIndelSize,
                                                    float covThreshold,
                                                    float percentIdentity,
                                                    float minKmerCov) {
        
        boolean corrected = false;
        int numKmers = kmers.size();
        
        int k = graph.getK();
        int expectedGapSize = k-1;
        
        ArrayList<Kmer> kmers2 = new ArrayList<>(numKmers + maxIndelSize);
        int numBadKmersSince = 0;
        Kmer kmer;
        for (int i=0; i<numKmers; ++i) {
            kmer = kmers.get(i);
            
            if (kmer.count >= covThreshold) {
                
                if (numBadKmersSince > 0) {
                    if (kmers2.isEmpty()) {
                        // use original sequence
                        for (int j=0; j<i; ++j) {
                            kmers2.add(kmers.get(j));
                        }
                    }
                    else if (numBadKmersSince == expectedGapSize) {
                        // a SNV bubble
                        Kmer leftEdgeKmer = kmers.get(i-numBadKmersSince);
                        Kmer rightEdgeKmer = kmers.get(i-1);
                        String prefix = leftEdgeKmer.toString().substring(0, k-1);
                        
                        ArrayList<Kmer> bestKmers = null;
                        float bestCov = Float.MIN_VALUE;
                        
                        ArrayList<Kmer> testKmers;
                        float[] m3;
                        for (Kmer var : rightEdgeKmer.getLeftVariants(k, graph.getMaxNumHash(), graph, minKmerCov)) {
                            testKmers = graph.getKmers(prefix + var.toString());
                            if (!testKmers.isEmpty()) {
                                m3 = getMinMedMaxKmerCoverage(testKmers);

                                float medCov = m3[1];
                                if (m3[0] >= minKmerCov && medCov > bestCov) {
                                    bestCov = medCov;
                                    bestKmers = testKmers;
                                }
                            }
                        }
                        
                        if (bestKmers != null & bestCov >= minKmerCov) {
                            // fill with best kmers
                            kmers2.addAll(bestKmers);
                            corrected = true;
                        }
                        else {
                            // fill with original kmers
                            for (int j=i-numBadKmersSince; j<i; ++j) {
                                kmers2.add(kmers.get(j));
                            }                            
                        }
                    }
                    else if (numBadKmersSince < expectedGapSize - maxIndelSize) {
                        ++numBadKmersSince;
                        continue;
                    }
                    else {
                        int maxLengthDifference = 0;
                        if (percentIdentity < 1.0f) {
                            maxLengthDifference = Math.max(1, (int) Math.ceil(numBadKmersSince * (1 - percentIdentity) / k)) * maxIndelSize;
                        }
                        
                        // find best candidate for left kmer
                        int kmers2Size = kmers2.size();
                        int bestLeftKmerIndex = kmers2Size-1;
                        int leftKmerStopIndex = Math.max(0, bestLeftKmerIndex-lookahead);
                        Kmer bestLeftKmer = kmers2.get(bestLeftKmerIndex);
                        for (int j=bestLeftKmerIndex-2; j>=leftKmerStopIndex; --j) {
                            Kmer tmpKmer = kmers2.get(j);
                            if (tmpKmer.count > bestLeftKmer.count) {
                                bestLeftKmer = tmpKmer;
                                bestLeftKmerIndex = j;
                            }
                        }
                        
                        // find best candidate for right kmer
                        int bestRightKmerIndex = i;
                        int rightKmerStopIndex = Math.min(bestRightKmerIndex+lookahead, kmers.size()-1);
                        Kmer bestRightKmer = kmer;
                        for (int j=bestRightKmerIndex+1; j<=rightKmerStopIndex; ++j) {
                            Kmer tmpKmer = kmers.get(j);
                            if (tmpKmer.count > bestRightKmer.count) {
                                bestRightKmer = tmpKmer;
                                bestRightKmerIndex = j;
                            }
                        }
                        
                        // find path between left and right kmers
                        ArrayDeque<Kmer> altPath = getMaxCoveragePath(graph, bestLeftKmer, bestRightKmer, numBadKmersSince + maxLengthDifference, lookahead, minKmerCov);
                        
                        if (altPath == null) {
                            // fill with original sequence
                            for (int j=i-numBadKmersSince; j<i; ++j) {
                                kmers2.add(kmers.get(j));
                            }
                        }
                        else {
                            float altPathMinCov = getMinimumKmerCoverage(altPath);
                            float oriPathMinCov = getMinimumKmerCoverage(kmers, i-numBadKmersSince, i);
                            
                            int altPathLen = altPath.size();

                            if (oriPathMinCov < altPathMinCov &&
                                    (oriPathMinCov < minKmerCov && altPathMinCov >= minKmerCov) ||
                                    (numBadKmersSince-maxLengthDifference <= altPathLen && altPathLen <= numBadKmersSince+maxLengthDifference && 
                                        (altPathLen <= k+maxIndelSize ||
                                            getPercentIdentity(graph.assemble(altPath), graph.assemble(kmers, i-numBadKmersSince, i)) >= percentIdentity))) {
                                
                                // backtrack to best left kmer
                                for (int j=kmers2.size()-1; j>bestLeftKmerIndex; --j) {
                                    kmers2.remove(j);
                                }
                                
                                // add all kmers from alternate path
                                kmers2.addAll(altPath);
                                
                                // set cursor to best right kmer
                                i = bestRightKmerIndex;
                                kmer = bestRightKmer;
                                
                                corrected = true;
                            }
                            else {
                                // fill with original sequence
                                for (int j=i-numBadKmersSince; j<i; ++j) {
                                    kmers2.add(kmers.get(j));
                                }
                            }
                        }
                    }

                    numBadKmersSince = 0;
                }
                
                kmers2.add(kmer);
            }
            else {
                ++numBadKmersSince;
            }
        }
        
        if (numBadKmersSince > 0 && numBadKmersSince < numKmers) {
            int i = numKmers-numBadKmersSince; // index of first bad kmer
            
            // use original sequence
            for (int j=i; j<numKmers; ++j) {
                kmers2.add(kmers.get(j));
            }
        }
        
        if (corrected) {
            return kmers2;
        }
        
        return null;
    }
    
    public static ArrayList<Kmer> correctErrorHelper(ArrayList<Kmer> kmers,
                                                    BloomFilterDeBruijnGraph graph, 
                                                    int lookahead,
                                                    int maxIndelSize,
                                                    float covThreshold,
                                                    float percentIdentity,
                                                    float minKmerCov) {
        
        
        boolean corrected = false;
        int numKmers = kmers.size();
        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        int expectedGapSize = k;
        
        ArrayList<Kmer> kmers2 = new ArrayList<>(numKmers + maxIndelSize);
        int numBadKmersSince = 0;
        Kmer kmer;
        for (int i=0; i<numKmers; ++i) {
            kmer = kmers.get(i);
            
            if (kmer.count >= covThreshold) {
                
                if (numBadKmersSince > 0) {
                    if (kmers2.isEmpty()) {
                        // check left edge
                        ArrayDeque<Kmer> leftVars = kmers.get(i-1).getLeftVariants(k, numHash, graph, minKmerCov);
                        if (leftVars.isEmpty()) {
                            // no branches found
                            for (int j=0; j<i; ++j) {
                                kmers2.add(kmers.get(j));
                            }
                        }
                        else if (numBadKmersSince >= lookahead) {
                            // calculate median cov of edge kmers
                            float[] tipCovs = new float[numBadKmersSince];
                            for (int j=0; j<i; ++j) {
                                tipCovs[j] = kmers.get(j).count;
                            }
                            float tipMedCov = getMedian(tipCovs);
                            
                            ArrayDeque<Kmer> greedyTipKmers = greedyExtendLeft(graph, kmer, lookahead, numBadKmersSince);
                            if (greedyTipKmers.size() == numBadKmersSince && getMedianKmerCoverage(greedyTipKmers) > tipMedCov) {
                                if (getPercentIdentity(graph.assemble(greedyTipKmers), graph.assemble(kmers, 0, i)) >= percentIdentity){
                                    corrected = true;
                                    kmers2.addAll(greedyTipKmers);
                                }
                                else if (!kmers.get(0).hasPredecessors(k, numHash, graph) && numBadKmersSince < k) {
                                    // this is blunt end in graph
                                    // do not add its kmers
                                    corrected = true;
                                }
                                else {
                                    // use original sequence
                                    for (int j=0; j<i; ++j) {
                                        kmers2.add(kmers.get(j));
                                    }
                                }
                            }
                            else {
                                // use original sequence
                                for (int j=0; j<i; ++j) {
                                    kmers2.add(kmers.get(j));
                                }
                            }
                        }
                        else {
                            corrected = true;
                        }
                    }
                    else if (numBadKmersSince == expectedGapSize) {
                        // a SNV bubble
                        Kmer leftEdgeKmer = kmers.get(i-numBadKmersSince);
                        Kmer rightEdgeKmer = kmers.get(i-1);
                        String left = leftEdgeKmer.toString();
                        String right = rightEdgeKmer.toString();
                        
                        ArrayList<Kmer> bestKmers = null;
                        float bestCov = Float.MIN_VALUE;
                        
                        ArrayList<Kmer> testKmers;
                        float[] m3;
                        for (char n : NUCLEOTIDES) {
                            testKmers = graph.getKmers(left + n + right);
                            if (!testKmers.isEmpty()) {
                                m3 = getMinMedMaxKmerCoverage(testKmers);

                                float medCov = m3[1];
                                if (m3[0] >= minKmerCov && medCov > bestCov) {
                                    bestCov = medCov;
                                    bestKmers = testKmers;
                                }
                            }
                        }
                        
                        if (bestKmers != null & bestCov >= minKmerCov) {
                            // fill with best kmers
                            kmers2.addAll(bestKmers);
                            corrected = true;
                        }
                        else {
                            // fill with original kmers
                            for (int j=i-numBadKmersSince; j<i; ++j) {
                                kmers2.add(kmers.get(j));
                            }                            
                        }
                    }
                    else {
                        ArrayDeque<Kmer> path = getMaxCoveragePath(graph, kmers2.get(kmers2.size()-1), kmer, numBadKmersSince + maxIndelSize, lookahead, minKmerCov);
                        if (path == null) {
                            // fill with original sequence
                            for (int j=i-numBadKmersSince; j<i; ++j) {
                                kmers2.add(kmers.get(j));
                            }
                        }
                        else {
                            int altPathLen = path.size();

                            if (numBadKmersSince-maxIndelSize <= altPathLen && 
                                    altPathLen <= numBadKmersSince+maxIndelSize && 
                                    (altPathLen <= k+maxIndelSize ||
                                        getPercentIdentity(graph.assemble(path), graph.assemble(kmers, i-numBadKmersSince, i)) >= percentIdentity)) {
                                
                                    kmers2.addAll(path);
                                    corrected = true;
                            }
                            else {
                                // fill with original sequence
                                for (int j=i-numBadKmersSince; j<i; ++j) {
                                    kmers2.add(kmers.get(j));
                                }
                            }
                        }
                    }

                    numBadKmersSince = 0;
                }

                kmers2.add(kmer);
            }
            else {
                ++numBadKmersSince;
            }
        }
        
        if (numBadKmersSince > 0 && numBadKmersSince < numKmers) {
            // check right edge
            int i = numKmers-numBadKmersSince; // index of first bad kmer
            ArrayDeque<Kmer> rightVars = kmers.get(i).getRightVariants(k, numHash, graph, minKmerCov);
            if (rightVars.isEmpty()) {
                for (int j=i; j<numKmers; ++j) {
                    kmers2.add(kmers.get(j));
                }
            }
            else if (numBadKmersSince >= lookahead) {
                // calculate median cov of edge kmers
                float[] tipCovs = new float[numBadKmersSince];
                for (int j=0; j<numBadKmersSince; ++j) {
                    tipCovs[j] = kmers.get(i+j).count;
                }
                float tipMedCov = getMedian(tipCovs);
                
                ArrayDeque<Kmer> greedyTipKmers = greedyExtendRight(graph, kmers.get(i-1), lookahead, numBadKmersSince);
                if (greedyTipKmers.size() == numBadKmersSince && getMedianKmerCoverage(greedyTipKmers) > tipMedCov) {
                    if (getPercentIdentity(graph.assemble(greedyTipKmers), graph.assemble(kmers, i, numKmers)) >= percentIdentity){
                        corrected = true;
                        kmers2.addAll(greedyTipKmers);
                    }
                    else if (!kmers.get(numKmers-1).hasSuccessors(k, numHash, graph) && numBadKmersSince < k) {
                        // this is blunt end in graph
                        // do not add its kmers
                        corrected = true;
                    }
                    else {
                        // use original sequence
                        for (int j=i; j<numKmers; ++j) {
                            kmers2.add(kmers.get(j));
                        }
                    }
                }
                else {
                    // use original sequence
                    for (int j=i; j<numKmers; ++j) {
                        kmers2.add(kmers.get(j));
                    }
                }
            }
            else {
                corrected = true;
            }
        }
        
        corrected = correctMismatches(kmers2, graph, covThreshold, minKmerCov) || corrected;
         
        if (corrected) {
            return kmers2;
        }
        
        // no changes
        return null;
    }
    
    public static boolean correctMismatches(ArrayList<Kmer> kmers,
                                            BloomFilterDeBruijnGraph graph,
                                            float covThreshold,
                                            float minKmerCov) {
        boolean corrected = false;
        int numKmers = kmers.size();
        int k = graph.getK();
        Kmer kmer;
        
        Kmer left, right;
        
        // scan in the forward direction
        for (int i=1; i<numKmers-k; ++i) {
            kmer = kmers.get(i);
            if (kmer.count < covThreshold) {
                left = kmers.get(i-1);
                if (left.count >= covThreshold) {
                    right = kmers.get(i+k);
                    if (right.count >= covThreshold) {
                        String tail = graph.getPrefix(right.toString());
                        ArrayList<Kmer> bestAlt = null;
                        float bestCov = getMedianKmerCoverage(kmers, i, i+k-1);
                        
                        for (String var : graph.getRightVariants(kmer.toString())) {
                            String alt = var + tail;
                            ArrayList<Kmer> altKmers = graph.getKmers(alt);
                            if (!altKmers.isEmpty()) {
                                float[] m = getMinMedMaxKmerCoverage(altKmers);
                                if (m[0] >= minKmerCov && m[1] > bestCov) {
                                    bestCov = m[1];
                                    bestAlt = altKmers;
                                }
                            }
                        }
                        
                        if (bestAlt != null) {
                            for (int j=0; j<k; ++j) {
                                kmers.set(i+j, bestAlt.get(j));
                            }
                            corrected = true;
                        }
                    }
                }
            }
        }
        
        // scan in the reverse direction
        for (int i=numKmers-2; i>=k; --i) {
            kmer = kmers.get(i);
            if (kmer.count < covThreshold) {
                right = kmers.get(i+1);
                if (right.count >= covThreshold) {
                    left = kmers.get(i-k);
                    if (left.count >= covThreshold) {
                        String head = graph.getSuffix(left.toString());
                        ArrayList<Kmer> bestAlt = null;
                        float bestCov = getMedianKmerCoverage(kmers, i-k+1, i);
                        
                        for (String var : graph.getLeftVariants(kmer.toString())) {
                            String alt = head + var;
                            ArrayList<Kmer> altKmers = graph.getKmers(alt);
                            if (!altKmers.isEmpty()) {
                                float[] m = getMinMedMaxKmerCoverage(altKmers);
                                if (m[0] >= minKmerCov && m[1] > bestCov) {
                                    bestCov = m[1];
                                    bestAlt = altKmers;
                                }
                            }
                        }
                        
                        if (bestAlt != null) {
                            for (int j=0; j<k; ++j) {
                                kmers.set(i-k+1+j, bestAlt.get(j));
                            }
                            corrected = true;
                        }
                    }
                }
            }
        }
        
        return corrected;
    }
            
    public static ArrayList<Kmer> correctErrorsSE(ArrayList<Kmer> kmers,
                                    BloomFilterDeBruijnGraph graph, 
                                    int lookahead,
                                    int maxIndelSize,
                                    float maxCovGradient, 
                                    float covFPR,
                                    float percentIdentity,
                                    float minKmerCov) {
        
        if (!kmers.isEmpty()) {
            int numKmers = kmers.size();
            
            // sort coverage in ascending order
            float[] covs = new float[numKmers];
            for (int i=0; i<numKmers; ++i) {
                covs[i] = kmers.get(i).count;
            }
            Arrays.sort(covs);

            boolean thresholdFound = false;

            int numFalsePositivesAllowed = (int) Math.round(numKmers * covFPR);
            int startIndex = numKmers - 1 - numFalsePositivesAllowed;

            float covThreshold = 0;

            if (startIndex >= 0) {
                covThreshold = covs[startIndex];
                float c = -1;
                for (int i=startIndex-1; i>=0; --i) {
                    c = covs[i];
                    if (covThreshold * maxCovGradient > c) {
                        thresholdFound = true;
                        break;
                    }
                    covThreshold = c;
                }
            }

            if (thresholdFound) {
                return correctErrorHelper(kmers,
                                            graph, 
                                            lookahead,
                                            maxIndelSize,
                                            covThreshold,
                                            percentIdentity,
                                            minKmerCov);
            }
        }

        return null;
    }
    
    public static ReadPair correctErrorsPE(ArrayList<Kmer> leftKmers, 
                                            ArrayList<Kmer> rightKmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead,
                                            int maxIndelSize, 
                                            float maxCovGradient, 
                                            float covFPR,
                                            int errorCorrectionIterations,
                                            float minCovThreshold,
                                            float percentIdentity,
                                            float minKmerCov) {
        
        boolean leftCorrected = false;
        boolean rightCorrected = false;
        
        for (int round=0; round<errorCorrectionIterations; ++round) {
        
            int numLeftKmers = leftKmers.size();
            int numRightKmers = rightKmers.size();

            int numFalsePositivesAllowed = (int) Math.round(Math.max(numLeftKmers, numRightKmers) * covFPR);

            // sort coverage of left kmers in ascending order
            float[] covs = new float[numLeftKmers];
            for (int i=0; i<numLeftKmers; ++i) {
                covs[i] = leftKmers.get(i).count;
            }
            Arrays.sort(covs);
            
//            float leftMedianCoverage = covs[numLeftKmers/2];

            // find cov threshold in left kmers
            boolean leftThresholdFound = false;
            int startIndex = numLeftKmers - 1;
            if (startIndex > numFalsePositivesAllowed) {
                startIndex -= numFalsePositivesAllowed;
            }
            float leftCovThreshold = covs[startIndex];
            float c;
            for (int i=startIndex-1; i>=0; --i) {
                c = covs[i];
                if (leftCovThreshold * maxCovGradient >= c) {
                    leftThresholdFound = true;
                    break;
                }
                leftCovThreshold = c;
            }

            // sort coverage of right kmers in ascending order
            covs = new float[numRightKmers];
            for (int i=0; i<numRightKmers; ++i) {
                covs[i] = rightKmers.get(i).count;
            }
            Arrays.sort(covs);

//            float rightMedianCov = covs[numRightKmers/2];
//            
//            if (Math.max(leftMedianCoverage, rightMedianCov) < medCovThreshold) {
//                break;
//            }
            
            // find cov threshold in right kmers
            boolean rightThresholdFound = false;
            startIndex = numRightKmers - 1;
            if (startIndex > numFalsePositivesAllowed) {
                startIndex -= numFalsePositivesAllowed;
            }
            float rightCovThreshold = covs[startIndex];
            for (int i=startIndex-1; i>=0; --i) {
                c = covs[i];
                if (rightCovThreshold * maxCovGradient >= c) {
                    rightThresholdFound = true;
                    break;
                }
                rightCovThreshold = c;
            }
            
            // set coverage threshold for both reads
            float covThreshold = -1;
            
            if (leftThresholdFound && rightThresholdFound) {
                covThreshold = Math.min(leftCovThreshold, rightCovThreshold);
            }
            else if (leftThresholdFound) {
                if (leftCovThreshold <= rightCovThreshold) {
                    covThreshold = leftCovThreshold;
                }
            }
            else if (rightThresholdFound) {
                if (rightCovThreshold <= leftCovThreshold) {
                    covThreshold = rightCovThreshold;
                }
            }
                
            if (covThreshold >= minCovThreshold) {
                // correct left read
                ArrayList<Kmer> leftKmers2 = correctErrorHelper(leftKmers,
                                                                graph, 
                                                                lookahead,
                                                                maxIndelSize,
                                                                covThreshold,
                                                                percentIdentity,
                                                                minKmerCov);

                if (leftKmers2 != null) {
                    leftKmers = leftKmers2;
                    leftCorrected = true;
                }

                // correct right read
                ArrayList<Kmer> rightKmers2 = correctErrorHelper(rightKmers,
                                                                graph, 
                                                                lookahead,
                                                                maxIndelSize,
                                                                covThreshold,
                                                                percentIdentity,
                                                                minKmerCov);

                if (rightKmers2 != null) {
                    rightKmers = rightKmers2;
                    rightCorrected = true;
                }

                if (leftKmers2 == null && rightKmers2 == null) {
                    break;
                }

            }
        }
        
        return new ReadPair(leftKmers, rightKmers, leftCorrected || rightCorrected);
    }
    
    public static ArrayDeque<int[]> breakWithReadPairedKmers(ArrayList<Kmer> kmers, BloomFilterDeBruijnGraph graph, int numPairsRequired, int rangeStart, int rangeEnd) {
        int interlockDistance = 0;
        
        ArrayDeque<int[]> segments = new ArrayDeque<>();
        
        int d = graph.getReadPairedKmerDistance();
        int lastIndex = rangeEnd - 1 - d;
        
        int start = -1;
        int end = -1;
        
        if (numPairsRequired == 1) {
            for (int i=rangeStart; i<=lastIndex; ++i) {
                if (graph.lookupReadKmerPair(kmers.get(i), kmers.get(i+d))) {
                    if (start < 0) {
                        start = i;
                    }

                    end = i+d;
                }
                else if (start >= 0 && i >= end-interlockDistance) {
                    segments.add(new int[]{start,end+1});
                    
                    start = -1;
                    end = -1;
                }
            }

            if (start >= 0) {
                segments.add(new int[]{start,end+1});
            }
        }
        else {
            int numPreviousKmerPairs = 0;
            for (int i=rangeStart; i<=lastIndex; ++i) {
                if (graph.lookupReadKmerPair(kmers.get(i), kmers.get(i+d))) {
                    if (++numPreviousKmerPairs >= numPairsRequired) {
                        if (start < 0) {
                            start = i-numPairsRequired+1;
                        }
                        
                        end = i+d;
                    }
                }
                else {
                    if (start >= 0 && i >= end-interlockDistance) {
                        segments.add(new int[]{start,end+1});

                        start = -1;
                        end = -1;
                    }
                    
                    numPreviousKmerPairs = 0;
                }
            }
            
            if (start >= 0) {
                segments.add(new int[]{start,end+1});
            }
        }
        
        return segments;
    }
    
    public static ArrayDeque<int[]> breakWithReadPairedKmers(ArrayList<Kmer> kmers, BloomFilterDeBruijnGraph graph, int numPairsRequired) {
        int interlockDistance = 0;
        
        ArrayDeque<int[]> segments = new ArrayDeque<>();
        
        int d = graph.getReadPairedKmerDistance();
        int lastIndex = kmers.size() - 1 - d;
        
        int start = -1;
        int end = -1;
        
        if (numPairsRequired == 1) {
            for (int i=0; i<=lastIndex; ++i) {
                if (graph.lookupReadKmerPair(kmers.get(i), kmers.get(i+d))) {
                    if (start < 0) {
                        start = i;
                    }

                    end = i+d;
                }
                else if (start >= 0 && i >= end-interlockDistance) {
                    segments.add(new int[]{start,end+1});
                    
                    start = -1;
                    end = -1;
                }
            }

            if (start >= 0) {
                segments.add(new int[]{start,end+1});
            }
        }
        else {
            int numPreviousKmerPairs = 0;
            for (int i=0; i<=lastIndex; ++i) {
                if (graph.lookupReadKmerPair(kmers.get(i), kmers.get(i+d))) {
                    if (++numPreviousKmerPairs >= numPairsRequired) {
                        if (start < 0) {
                            start = i-numPairsRequired+1;
                        }
                        
                        end = i+d;
                    }
                }
                else {
                    if (start >= 0 && i >= end-interlockDistance) {
                        segments.add(new int[]{start,end+1});

                        start = -1;
                        end = -1;
                    }
                    
                    numPreviousKmerPairs = 0;
                }
            }
            
            if (start >= 0) {
                segments.add(new int[]{start,end+1});
            }
        }
        
        return segments;
    }
    
    public static ArrayDeque<int[]> breakWithFragPairedKmers(ArrayList<Kmer> kmers, BloomFilterDeBruijnGraph graph, int numPairsRequired) {
        int interlockDistance = 0;
        
        ArrayDeque<int[]> segments = new ArrayDeque<>();
        
        int d = graph.getFragPairedKmerDistance();
        int lastIndex = kmers.size() - 1 - d;
        
        int start = -1;
        int end = -1;
        
        if (numPairsRequired == 1) {
            for (int i=0; i<=lastIndex; ++i) {
                if (graph.lookupKmerPair(kmers.get(i), kmers.get(i+d))) {
                    if (start < 0) {
                        start = i;
                    }

                    end = i+d;
                }
                else if (start >= 0 && i >= end-interlockDistance) {
                    segments.add(new int[]{start, end+1});

                    start = -1;
                    end = -1;
                }
            }

            if (start >= 0) {
                segments.add(new int[]{start, end+1});
            }
        }
        else {
            int numPreviousKmerPairs = 0;
            for (int i=0; i<=lastIndex; ++i) {
                if (graph.lookupKmerPair(kmers.get(i), kmers.get(i+d))) {
                    if (++numPreviousKmerPairs >= numPairsRequired) {
                        if (start < 0) {
                            start = i-numPairsRequired+1;
                        }
                        
                        end = i+d;
                    }
                }
                else {
                    if (start >= 0 && i >= end-interlockDistance) {
                        segments.add(new int[]{start, end+1});

                        start = -1;
                        end = -1;
                    }
                    
                    numPreviousKmerPairs = 0;
                }
            }
            
            if (start >= 0) {
                segments.add(new int[]{start, end+1});
            }
        }
        
        return segments;
    }
    
    public static ArrayDeque<int[]> breakWithFragPairedKmers(ArrayList<Kmer> kmers, BloomFilterDeBruijnGraph graph) {
        int interlockDistance = 0;
        
        ArrayDeque<int[]> segments = new ArrayDeque<>();
        
        int d = graph.getFragPairedKmerDistance();
        int lastIndex = kmers.size() - 1 - d;
        
        int start = -1;
        int end = -1;
        
        for (int i=0; i<=lastIndex; ++i) {
            if (graph.lookupKmerPair(kmers.get(i), kmers.get(i+d))) {
                if (start < 0) {
                    start = i;
                }
                
                end = i+d;
            }
            else if (start >= 0 && i >= end-interlockDistance) {
                segments.add(new int[]{start, end+1});
                
                start = -1;
                end = -1;
            }
        }
        
        if (start >= 0) {
            segments.add(new int[]{start, end+1});
        }
        
        return segments;
    }
    
//    public static String correctErrors(String seq, BloomFilterDeBruijnGraph graph, int lookahead, int errorsAllowed) {
//        int numCorrected = 0;
//        int seqLen = seq.length();
//        int k = graph.getK();
//        
//        if (seqLen < k) {
//            // no correction
//            return seq;
//        }
//        
//        int numKmers = seqLen-k+1;
//        StringBuilder sb = new StringBuilder(seq);
//        
//        float bestCov, cov;
//        String kmer, guide;
//        ArrayDeque<String> variants;
//        
//        // correct from start
//        for (int i=0; i<numKmers; ++i) {
//            int j = i+k;
//            kmer = sb.substring(i, j);
//            variants = graph.getRightVariants(kmer);
//            if (!variants.isEmpty()) {
//                guide = sb.substring(j, Math.min(j+lookahead, seqLen));
//                int guideLen = guide.length();
//                bestCov = 0;
//                
//                if (graph.contains(kmer)) {
//                    bestCov = rightGuidedMedianCoverage(graph, kmer, guide, guideLen);
//                }
//                
//                // test for mismatch
//                boolean correctedMismatch = false;
//                char bestBase = 'N';
//                for (String v : variants) {
//                    cov = rightGuidedMedianCoverage(graph, v, guide, guideLen);
//                                        
//                    if (cov > bestCov) {
//                        bestCov = cov;
//                        bestBase = v.charAt(k-1);
//                        correctedMismatch = true;
//                    }
//                }
//                
//                // test for insertion in the sequence, ie. last base of kmer is the inserted base
//                boolean correctedInsertion = false;
//                if (j < seqLen-2) {
//                    String a = graph.getPrefix(kmer) + sb.charAt(j);
//                    if (graph.contains(a)) {
//                        String guideIns = sb.substring(j+1, Math.min(j+1+lookahead, seqLen));
//
//                        cov = rightGuidedMedianCoverage(graph, a, guideIns, guideIns.length());
//                        if (cov > bestCov) {
//                            correctedMismatch = false;
//                            correctedInsertion = true;
//
//                            bestCov = cov;
//                        }
//                    }
//                }
//                
//                // test for deletion in the sequence
//                boolean correctedDeletion = false;
//                String guideDel = sb.substring(j-1, Math.min(j-1+lookahead, seqLen));
//                int guideDelLen = guideDel.length();
//                char bestInsBase = 'N';
//                for (String v : variants) {
//                    cov = rightGuidedMedianCoverage(graph, v, guideDel, guideDelLen);
//                                        
//                    if (cov > bestCov) {
//                        correctedMismatch = false;
//                        correctedInsertion = false;
//                        correctedDeletion = true;
//                        
//                        bestCov = cov;
//                        bestInsBase = v.charAt(k-1);
//                    }
//                }
//                
//
//                if (correctedMismatch) {
//                    if (++numCorrected > errorsAllowed) {
//                        // too many errors; do not apply corrections
//                        return seq;
//                    }
//
//                    // replace the mismatch base
//                    sb.setCharAt(j-1, bestBase);
//                }
//                else if (correctedInsertion) {
//                    if (++numCorrected > errorsAllowed) {
//                        // too many errors; do not apply corrections
//                        return seq;
//                    }
//
//                    // remove the inserted base
//                    sb.deleteCharAt(j-1);
//                    --seqLen;
//                    --numKmers;
//                }
//                else if (correctedDeletion) {
//                    if (++numCorrected > errorsAllowed) {
//                        // too many errors; do not apply corrections
//                        return seq;
//                    }
//
//                    // insert the deleted base
//                    sb.insert(j-1, bestInsBase);
//                    ++seqLen;
//                    ++numKmers;
//                }
//            }
//        }
//        
//        // correct from end
//        for (int i=seqLen-k; i>=0; --i) {
//            kmer = sb.substring(i, i+k);
//            variants = graph.getLeftVariants(kmer);
//            if (!variants.isEmpty()) {
//                guide = sb.substring(Math.max(0, i-lookahead), i);
//                int guideLen = guide.length();
//                bestCov = 0;
//                
//                if (graph.contains(kmer)) {
//                    bestCov = leftGuidedMedianCoverage(graph, kmer, guide, guideLen);
//                }
//                                
//                // test for mismatch
//                boolean correctedMismatch = false;
//                char bestBase = 'N';
//                for (String v : variants) {
//                    cov = leftGuidedMedianCoverage(graph, v, guide, guideLen);
//                                        
//                    if (cov > bestCov) {
//                        bestCov = cov;
//                        bestBase = v.charAt(0);
//                        correctedMismatch = true;
//                    }
//                }
//                
//                // test for insertion in the sequence; ie. first base of kmer is the inserted base
//                boolean correctedInsertion = false;
//                if (i > 0) {
//                    String a = sb.charAt(i-1) + graph.getSuffix(kmer);
//                    if (graph.contains(a)) {
//                        String guideIns = sb.substring(Math.max(0, i-1-lookahead), i-1);
//
//                        cov = leftGuidedMedianCoverage(graph, a, guideIns, guideIns.length());
//                        if (cov > bestCov) {
//                            correctedMismatch = false;
//                            correctedInsertion = true;
//
//                            bestCov = cov;
//                        }
//                    }
//                }
//                
//                // test for deletion in the sequence
//                boolean correctedDeletion = false;
//                String guideDel = sb.substring(Math.max(0, i+1-lookahead), i+1);
//                int guideDelLen = guideDel.length();
//                char bestInsBase = 'N';
//                for (String v : variants) {
//                    cov = leftGuidedMedianCoverage(graph, v, guideDel, guideDelLen);
//                                        
//                    if (cov > bestCov) {
//                        correctedMismatch = false;
//                        correctedInsertion = false;
//                        correctedDeletion = true;
//                        
//                        bestCov = cov;
//                        bestInsBase = v.charAt(0);
//                    }
//                }
//                
//                if (correctedMismatch) {
//                    if (++numCorrected > errorsAllowed) {
//                        // too many errors; do not apply corrections
//                        return seq;
//                    }
//
//                    // replace the mismatch base
//                    sb.setCharAt(i, bestBase);
//                }
//                else if (correctedInsertion) {
//                    if (++numCorrected > errorsAllowed) {
//                        // too many errors; do not apply corrections
//                        return seq;
//                    }
//
//                    // remove the inserted base
//                    sb.deleteCharAt(i);
//                    --seqLen;
//                    --numKmers;
//                }
//                else if (correctedDeletion) {
//                    if (++numCorrected > errorsAllowed) {
//                        // too many errors; do not apply corrections
//                        return seq;
//                    }
//
//                    // insert the deleted base
//                    sb.insert(i, bestInsBase);
//                    ++seqLen;
//                    ++numKmers;
//                }
//            }
//        }
//        
//        if (numCorrected == 0) {
//            return seq;
//        }
//        
//        String seq2 = sb.toString();
//        if (graph.isValidSeq(seq2)) {
//            return seq2;
//        }
//        
//        return seq;
//    }
    
//    public static float[] coverageGradients(String seq, BloomFilterDeBruijnGraph graph, int lookahead) {
//        float[] counts = graph.getCounts(seq);
//        int numCounts = counts.length;
//        
//        ArrayDeque<Float> window = new ArrayDeque<>();
//        for (int i=0; i<lookahead; ++i) {
//            window.addLast(counts[i]);
//        }        
//
//        int numMedCounts = numCounts-lookahead+1;
//        float[] medCounts = new float[numMedCounts];
//        medCounts[0] = getMedian(window);
//        int m = 0;
//        for (int i=lookahead; i<numCounts; ++i) {
//            window.removeFirst();
//            window.addLast(counts[i]);
//            medCounts[++m] = getMedian(window);
//        }
//        
//        int numGradients = numCounts-(2*lookahead)+1;
//        float[] gradients = new float[numGradients];
//        for (int i=0; i<numGradients; ++i) {
//            float r = medCounts[i]/medCounts[i+lookahead];
//            if (r > 1) {
//                gradients[i] = 1/r;
//            }
//            else {
//                gradients[i] = r;
//            }
//        }
//        
//        return gradients;
//    }

    public static String assembleString(ArrayList<String> kmers) {
        
        String first = kmers.get(0);
        int k = first.length();
        int lastIndex = k - 1;
        
        StringBuilder sb = new StringBuilder(k + kmers.size() - 1);
        sb.append(first.substring(0, lastIndex));
        
        for (String kmer : kmers) {
            sb.append(kmer.charAt(lastIndex));
        }
        
        return sb.toString();
    }
    
    public static String assembleKmers(String[] kmers, int k) {
        String first = kmers[0];
        int lastIndex = k - 1;
        
        StringBuilder sb = new StringBuilder(k + kmers.length - 1);
        sb.append(first.substring(0, lastIndex));
        
        for (String kmer : kmers) {
            sb.append(kmer.charAt(lastIndex));
        }
        
        return sb.toString();
    }
        
    public static ArrayList<Kmer> greedyExtend(Kmer seed, BloomFilterDeBruijnGraph graph, int lookahead) {
        /**@TODO store smallest strand kmers for non-strand specific sequences */
                
        HashSet<String> pathKmerStr = new HashSet<>(1000);
        pathKmerStr.add(seed.toString());
        
        ArrayList<Kmer> rightPath = new ArrayList<>(1000);
        
        /* extend on right side */
        Kmer best = seed;
        while (true) {
            best = greedyExtendRightOnce(graph, best, lookahead);
            if (best != null) {
                String seq = best.toString();
                if (pathKmerStr.contains(seq)) {
                    break;
                }
                pathKmerStr.add(seq);
                rightPath.add(best);
            }
            else {
                break;
            }
        }
        
        ArrayList<Kmer> leftPath = new ArrayList<>(100);
        
        /* extend on left side */
        best = seed;
        while (true) {
            best = greedyExtendLeftOnce(graph, best, lookahead);
            if (best != null) {
                String seq = best.toString();
                if (pathKmerStr.contains(seq)) {
                    break;
                }
                pathKmerStr.add(seq);
                leftPath.add(best);
            }
            else {
                break;
            }
        }
        
        Collections.reverse(leftPath);
        leftPath.add(seed);
        leftPath.addAll(rightPath);
        
        return leftPath;
    }
    
    public static Kmer findMaxCoverageWindowKmer(ArrayList<Kmer> path, BloomFilterDeBruijnGraph graph, int windowSize) {
        int pathLen = path.size();
        if (pathLen <= windowSize) {
            return path.get(pathLen/2);
        }
                
        ArrayDeque<Float> window = new ArrayDeque<>();
        for (int i=0; i<windowSize; ++i) {
            window.addFirst(path.get(i).count);
        }
        
        int start = 0;
        int end = windowSize;
        float maxCov = getMedian(window);
        
        float cov;
        for (int i=windowSize; i<pathLen; ++i) {
            window.removeLast();
            window.addFirst(path.get(i).count);
            cov = getMedian(window);
            if (cov > maxCov) {
                maxCov = cov;
                end = i+1;
                start = end-windowSize;
            }
        }
        
        return path.get((end+start)/2);
    }
    
    public static ArrayList<Kmer> findBackbonePath(Kmer seed, BloomFilterDeBruijnGraph graph, int lookahead, int windowSize, int maxIteration) {
        Kmer best = seed;
        ArrayList<Kmer> path = greedyExtend(best, graph, lookahead);
        boolean randomSeed = false;
        
        for (int i=1; i<maxIteration; ++i) {
            if (randomSeed) {
                best = path.get((int) (Math.random() * (path.size()-1)));
                randomSeed = false;
            }
            else {
                best = findMaxCoverageWindowKmer(path, graph, windowSize);
                randomSeed = true;
            }
            path = greedyExtend(best, graph, lookahead);
        }
        
        return path;
    }
    
    public static String getBestSegment(ArrayList<String> segments, BloomFilterDeBruijnGraph graph) {
        int numSeqs = segments.size();
        switch (numSeqs) {
            case 0:
                return "";
            case 1:
                return segments.get(0);
            default:
                int k = graph.getK();
                
                String best = "";
                float bestCov = 0;
                
                NTHashIterator itr = graph.getHashIterator();
                long[] hVals = itr.hVals;
                
                for (String seg : segments) {
                    if (seg.length() >= k) {
                        itr.start(seg);
                        
                        float min = Float.MAX_VALUE;
                        while (itr.hasNext()) {
                            itr.next();
                            float c = graph.getCount(hVals);
                            if (c < min) {
                                min = c;
                            }
                        }
                        
                        if (min > bestCov) {
                            best = seg;
                            bestCov = min;
                        }
                    }
                }
                
                return best;
        }        
    }
    
    public static String connect(ArrayList<String> segments, BloomFilterDeBruijnGraph graph, int lookahead) {
        int numSeqs = segments.size();
        switch (numSeqs) {
            case 0:
                return "";
            case 1:
                return segments.get(0);
            default:
                int k = graph.getK();
                String last = segments.get(0);
                String longest = last;
    
                for (int i=1; i<numSeqs; i+=2) {             
                    String current = segments.get(i+1);     
                    
                    String connected = connect(last, current, graph, segments.get(i).length()+k, lookahead);
                    int connectedLength = connected.length();
                    
                    if (connectedLength > 0) {
                        last = connected;
                        
                        if (connectedLength > longest.length()) {
                            longest = connected;
                        }
                    }
                    else {
                        last = current;
                        
                        if (current.length() > longest.length()) {
                            longest = current;
                        }
                    }
                }
                
                return longest;
        }
    }
    
    public static String connect(String left, String right, BloomFilterDeBruijnGraph graph, int bound, int lookahead) {
        int k = graph.getK();
        
        ArrayDeque<Kmer> pathKmers = getMaxCoveragePath(graph, graph.getKmer(getLastKmer(left, k)), graph.getKmer(getFirstKmer(right, k)), bound, lookahead, 1);
        
        if (pathKmers == null || pathKmers.isEmpty() || pathKmers.size() > bound) {
            return "";
        }
        
        String leftWing, rightWing;
        int leftReadLength = left.length();
        if (leftReadLength == k) {
            // first base only
            leftWing = left.substring(0, 1);
        }
        else {
            leftWing = left.substring(0, leftReadLength-k+1);
        }
        
        rightWing = right.substring(k-1);
        
        return leftWing + graph.assemble(pathKmers) + rightWing;
    }

    public static ArrayList<Kmer> overlap(ArrayList<Kmer> leftKmers, ArrayList<Kmer> rightKmers, BloomFilterDeBruijnGraph graph, int minOverlap, float minKmerCov) {        
        String left = graph.assemble(leftKmers);
        String right = graph.assemble(rightKmers);
        String overlapped = overlapMaximally(left, right, minOverlap);
        
        if (overlapped == null) {
            // detect dovetail overlap
            
            // increase min overlap threshold
            minOverlap = Math.max(minOverlap, Math.min(left.length(), right.length())*3/4);
            
            overlapped = overlapMaximally(right, left, minOverlap);
            
            if (overlapped != null) {
                //System.out.println(">left\n" + left + "\n>right\n" + right + "\n>overlapped\n" + overlapped);
                
                // swap left and right reads
                ArrayList<Kmer> tmpKmers = leftKmers;
                leftKmers = rightKmers;
                rightKmers = tmpKmers;
                
                String tmp = left;
                left = right;
                right = tmp;
            }
        }
        
        if (overlapped != null) {
            int overlappedSeqLength = overlapped.length();
            int k = graph.getK();
                        
            int leftLen = left.length();
            int rightLen = right.length();
            
            if (overlappedSeqLength <= leftLen + rightLen - k) {
                // The overlap is larger than or equal to k
                
                if (overlappedSeqLength == Math.max(leftLen, rightLen)) {
                    if (leftLen >= rightLen) {
                        // left read contains right read
                        return leftKmers;
                    }
                    else {
                        // right read contains left read
                        return rightKmers;
                    }
                }
                else {
                    boolean hasComplexKmer = false;

                    int end = rightLen - (overlappedSeqLength - leftLen) - k + 1;

                    for (int i=0; i<end; ++i) {
                        if (!isHomopolymer(rightKmers.get(i).bytes)) {
                            // Require at least one complex kmers in the overlap
                            hasComplexKmer = true;
                            break;
                        }
                    }
                    
                    if (!hasComplexKmer) {
                        return null;
                    }

                    ArrayList<Kmer> overlappedKmers = new ArrayList<>(overlappedSeqLength - k + 1); //graph.getKmers(overlapped);
                    overlappedKmers.addAll(leftKmers);

                    // add remaining right kmers
                    int numRightKmers = rightKmers.size();
                    for (int i=end; i<numRightKmers; ++i) {
                        overlappedKmers.add(rightKmers.get(i));
                    }
                    
                    return overlappedKmers;
                }
            }
            else {
                // The overlap is smaller than k
                boolean spanningKmersNotValid = false;
                boolean hasComplexKmer = false;
                boolean addReadKmers = false;
                int start = leftLen - k + 1;
                int end = overlappedSeqLength - (rightLen - k + 1);

                ArrayList<Kmer> spanningKmers = graph.getKmers(overlapped, start, end);
                for (Kmer kmer : spanningKmers) {
                    if (kmer.count < minKmerCov) {
                        // the overlap is not a valid path in the graph
                        spanningKmersNotValid = true;
                        break;
                    }
                    
                    if (!hasComplexKmer && !isHomopolymer(kmer.bytes)) {
                        // Require at least one complex kmers in the overlap
                        hasComplexKmer = true;
                    }
                }
                
                if (spanningKmersNotValid) {
                    // check whether kmers at edges of left and right reads have count = 1
                    int numBasesOverlapped = leftLen + rightLen - overlappedSeqLength;
                    boolean singletonInLeftRead = false;
                    boolean singeltonInRightRead = false;
                    int tmp = Math.min(numBasesOverlapped, rightKmers.size());
                    for (int i=0; i<tmp; ++i) {
                        if (rightKmers.get(i).count == 1) {
                            singeltonInRightRead = true;
                            break;
                        }
                    }
                    
                    if (singeltonInRightRead) {
                        int numLeftKmers = leftKmers.size();
                        for (int i=Math.max(0, numLeftKmers-numBasesOverlapped); i<numLeftKmers; ++i) {
                            if (leftKmers.get(i).count == 1) {
                                singletonInLeftRead = true;
                                break;
                            }
                        }
                        
                        if (singletonInLeftRead && !isRepeat(right.substring(0, numBasesOverlapped))) {
                            // add the missing kmers
                            for (Kmer kmer : spanningKmers) {
                                if (kmer.count == 0) {
                                    graph.addDbgOnly(kmer.getHash());
                                    kmer.count = 1;
                                }
                                
                                if (!hasComplexKmer && !isHomopolymer(kmer.bytes)) {
                                    // Require at least one complex kmers in the overlap
                                    hasComplexKmer = true;
                                }
                            }
                            
                            // add the missing read paired kmers
                            addReadKmers = true;
                        } 
                        else {
                            return null;
                        }
                    }
                    else {
                        return null;
                    }
                }
                
                if (!spanningKmers.isEmpty() && !hasComplexKmer) {
                    return null;
                }
                
                ArrayList<Kmer> overlappedKmers = new ArrayList<>(overlappedSeqLength - k + 1); //graph.getKmers(overlapped);
                overlappedKmers.addAll(leftKmers);
                overlappedKmers.addAll(spanningKmers);
                overlappedKmers.addAll(rightKmers);
                
                if (addReadKmers) {
                    correctMismatches(overlappedKmers, graph, 2, minKmerCov);
                    graph.addReadPairedKmers(overlappedKmers);
                }
                
                return overlappedKmers;
            }
        }
        
        return null;
    }
    
    public static ArrayList<Kmer> overlapAndConnect(ArrayList<Kmer> leftKmers, 
                                                    ArrayList<Kmer> rightKmers, 
                                                    BloomFilterDeBruijnGraph graph,
                                                    int bound, 
                                                    int lookahead,
                                                    int minOverlap,
                                                    float maxCovGradient,
                                                    int maxTipLen,
                                                    int maxIndelLen,
                                                    float minPercentIdentity,
                                                    float minKmerCov) {

        // 1. Attempt simple overlap
        ArrayList<Kmer> fragmentKmers = overlap(leftKmers, rightKmers, graph, minOverlap, minKmerCov);
        
        if (fragmentKmers == null) {

            fragmentKmers = join(graph, leftKmers, rightKmers, bound, lookahead, maxCovGradient,
                                    maxTipLen, maxIndelLen, minPercentIdentity, minKmerCov);
//            fragmentKmers = getSimilarCoveragePath(graph, leftKmers, rightKmers, bound, lookahead, maxCovGradient, false);
//            ArrayDeque<Kmer> connectedPath = getMaxCoveragePath(graph, leftLastKmer, rightFirstKmer, bound, lookahead);

        }
        
        return fragmentKmers;
    }
    
    public static ArrayList<Kmer> connect(final ArrayList<Kmer> leftKmers, 
                                            final ArrayList<Kmer> rightKmers, 
                                            final BloomFilterDeBruijnGraph graph,
                                            final int bound, 
                                            final int lookahead,
                                            final int minOverlap,
                                            final float maxCovGradient,
                                            final boolean rescueSearch) {
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        
        int numLeftKmers = leftKmers.size();
        int numRightKmers = rightKmers.size();
        
        HashSet<Kmer> leftKmersSet = new HashSet<>(leftKmers);
        HashSet<Kmer> rightKmersSet = new HashSet<>(rightKmers);
        
        // extend to RIGHT from left kmers
        
        ArrayList<Kmer> leftPathKmers = new ArrayList<>();
        HashSet<Kmer> leftPathKmersSet = new HashSet<>();
        
        ArrayDeque<Kmer> neighbors = leftKmers.get(leftKmers.size()-1).getSuccessors(k, numHash, graph);
        int depth = 0;
        int rightStart = -1;
        
        while (!neighbors.isEmpty() && depth <= bound) {
            Kmer cursor = null;
            
            if (neighbors.size() == 1) {
                cursor = neighbors.pop();
            }
            else {
                Kmer best = null;
                Kmer second = null;
                
                for (Kmer kmer : neighbors) {
                    if (best == null) {
                        best = kmer;
                    }
                    else if (second == null) {
                        if (best.count < kmer.count) {
                            second = best;
                            best = kmer;
                        }
                        else {
                            second = kmer;
                        }
                    }
                    else {
                        if (best.count < kmer.count) {
                            second = best;
                            best = kmer;
                        }
                        else if (second.count < kmer.count) {
                            second = kmer;
                        }                        
                    }
                    
                    if (rightKmersSet.contains(kmer)) {
                        for (rightStart=0; rightStart<numRightKmers; ++rightStart) {
                            if (kmer.equals(rightKmers.get(rightStart))) {
                                break;
                            }
                        }

                        if (rightStart == 0) {
                            ArrayList<Kmer> path = new ArrayList<>(leftKmers.size() + leftPathKmers.size() + rightKmers.size());
                            path.addAll(leftKmers);
                            path.addAll(leftPathKmers);
                            path.addAll(rightKmers);
                            return path;
                        }
                        else {
                            // Need to evaluate the dangling kmers later
                            break;
                        }
                    }
                }
                
                if (best.count * maxCovGradient >= second.count) {
                    cursor = best;
                }
            }
            
            if (cursor == null) {
                break;
            }
            
            if (rightKmersSet.contains(cursor)) {
                for (rightStart=0; rightStart<numRightKmers; ++rightStart) {
                    if (cursor.equals(rightKmers.get(rightStart))) {
                        break;
                    }
                }
                
                if (rightStart == 0) {
                    ArrayList<Kmer> path = new ArrayList<>(leftKmers.size() + leftPathKmers.size() + rightKmers.size());
                    path.addAll(leftKmers);
                    path.addAll(leftPathKmers);
                    path.addAll(rightKmers);
                    return path;
                }
                else {
                    // Need to evaluate the dangling kmers later
                    break;
                }
            }
            
            if (leftPathKmersSet.contains(cursor) || leftKmersSet.contains(cursor)) {
                // a loop
                break;
            }
            
            leftPathKmers.add(cursor);
            leftPathKmersSet.add(cursor);
            neighbors = cursor.getSuccessors(k, numHash, graph);
            ++depth;
        }
        
        // extend to LEFT from right kmers
        
        ArrayList<Kmer> rightPathKmers = new ArrayList<>();
        HashSet<Kmer> rightPathKmersSet = new HashSet<>();
        int leftStart = -1;
        
        depth = 0;
        neighbors = rightKmers.get(0).getPredecessors(k, numHash, graph);
        while (!neighbors.isEmpty() && depth <= bound) {
            Kmer cursor = null;
            
            if (neighbors.size() == 1) {
                cursor = neighbors.pop();
            }
            else {
                Kmer best = null;
                Kmer second = null;
                
                for (Kmer kmer : neighbors) {
                    if (best == null) {
                        best = kmer;
                    }
                    else if (second == null) {
                        if (best.count < kmer.count) {
                            second = best;
                            best = kmer;
                        }
                        else {
                            second = kmer;
                        }
                    }
                    else {
                        if (best.count < kmer.count) {
                            second = best;
                            best = kmer;
                        }
                        else if (second.count < kmer.count) {
                            second = kmer;
                        }                        
                    }
                }
                
                if (best.count * maxCovGradient >= second.count) {
                    cursor = best;
                }
            }
            
            if (cursor == null) {
                break;
            }
            
            if (leftKmersSet.contains(cursor)) {
                for (leftStart=0; leftStart<numLeftKmers; ++leftStart) {
                    if (cursor.equals(leftKmers.get(leftStart))) {
                        break;
                    }
                }
                
                if (leftStart == numLeftKmers-1) {
                    ArrayList<Kmer> path = new ArrayList<>(leftKmers.size() + rightPathKmers.size() + rightKmers.size());
                    path.addAll(leftKmers);
                    for (int i=rightPathKmers.size()-1; i>=0; --i) {
                        path.add(rightPathKmers.get(i));
                    }
                    path.addAll(rightKmers);
                    
                    return path;
                }
                else {
                    // Need to evaluate the dangling kmers later
                    break;
                }
            }
            else if (leftPathKmersSet.contains(cursor)) {
                // left path and right path intersect
                
                ArrayList<Kmer> path = new ArrayList<>(leftKmers.size() + leftPathKmers.size() + rightPathKmers.size() + rightKmers.size());
                path.addAll(leftKmers);
                for (Kmer kmer : leftPathKmers) {
                    if (kmer.equals(cursor)) {
                        break;
                    }
                    path.add(kmer);
                }
                path.add(cursor);
                for (int i=rightPathKmers.size()-1; i>=0; --i) {
                    path.add(rightPathKmers.get(i));
                }
                path.addAll(rightKmers);
                
                return path;
            }
            
            if (rightPathKmersSet.contains(cursor) || rightKmersSet.contains(cursor)) {
                // a loop
                break;
            }
            
            rightPathKmers.add(cursor);
            rightPathKmersSet.add(cursor);
            neighbors = cursor.getPredecessors(k, numHash, graph);
            ++depth;
        }
        
        if (leftStart > 0) {
            // TODO
        }
        
        if (rightStart > 0) {
            // TODO
        }
        
        return null;
    }

//    public static String extendWithPairedKmersNonGreedily(String fragment, BloomFilterDeBruijnGraph graph, int lookahead) {
//        final int distance = graph.getPairedKmerDistance();
//        final int k = graph.getK();
//        
//        // transcript kmers list
//        final ArrayList<String> kmers = new ArrayList<>(2*(fragment.length()-k+1));
//        kmerizeToCollection(fragment, k, kmers);
//        
//        // transcript kmers set
//        final HashSet<String> transcriptKmers = new HashSet<>(kmers);
//        
//        // kmer pairs used in the extension of this transcript
//        final HashSet<String> usedPairs = new HashSet<>();
//        
//        /*
//        float fragMinCov = Float.POSITIVE_INFINITY;
//        for (String kmer : kmers) {
//            float c = graph.getCount(kmer);
//            if (c < fragMinCov) {
//                fragMinCov = c;
//            }
//        }
//        //System.out.println(fragMinCov);
//        */
//        
//        /** extend right*/
//        String best = kmers.get(kmers.size()-1);
//        LinkedList<String> neighbors = graph.getSuccessors(best);
//        while (!neighbors.isEmpty()) {
//            best = null;
//            String partner = null;
//            
//            if (neighbors.size() == 1) {
//                // 1 successor
//                best = neighbors.peek();
//            }
//            else {
//                // >1 successors
//                
//                LinkedList<String> fragmentNeighbors = new LinkedList<>();
//                for (String n : neighbors) {
//                    if (graph.lookupFragmentKmer(n)) {
//                        fragmentNeighbors.add(n);
//                    }
//                }
//                
//                if (fragmentNeighbors.isEmpty()) {
//                    // 0 fragment successors
//                    break;
//                }
//                else if (fragmentNeighbors.size() == 1) {
//                    // 1 fragment successor
//                    best = fragmentNeighbors.peek();
//                }
//                else {
//                    // >1 fragment successors
//                    int numKmers = kmers.size();
//                    
//                    if (numKmers >= distance) {
//                        LinkedList<String> pairedNeighbors = new LinkedList<>();
//                        
//                        partner = kmers.get(numKmers-distance);
//                        for (String n : fragmentNeighbors) {
//                            if (graph.lookupPairedKmers(partner, n)) {
//                                pairedNeighbors.add(n);
//                            }
//                        }
//                        
//                        if (pairedNeighbors.size() == 1) {
//                            best = pairedNeighbors.peek();
//                        }
//                    }
//                }
//            }
//            
//            if (best == null) {
//                break;
//            }
//            else {
//                int numKmers = kmers.size();
//                
//                if (numKmers >= distance) {
//                    if (partner == null) {
//                        partner = kmers.get(numKmers-distance);
//                    }
//                    
//                    if (transcriptKmers.contains(best)) {
//                        if (partner != null) {
//                            String pairedKmersStr = partner + best;
//                            if (!graph.lookupPairedKmers(partner, best) || usedPairs.contains(pairedKmersStr)) {
//                                break;
//                            }
//                            else {
//                                usedPairs.add(pairedKmersStr);
//                            }
//                        }
//                    }
//                }
//            }
//            
//            kmers.add(best);
//            transcriptKmers.add(best);
//            neighbors = graph.getSuccessors(best);
//        }
//        
//        Collections.reverse(kmers);
//        
//        /** extend left*/
//        best = kmers.get(kmers.size()-1);
//        neighbors = graph.getPredecessors(best);
//        while (!neighbors.isEmpty()) {
//            best = null;
//            String partner = null;
//            
//            if (neighbors.size() == 1) {
//                best = neighbors.peek();
//            }
//            else {
//                // >1 neighbors
//                
//                LinkedList<String> fragmentNeighbors = new LinkedList<>();
//                for (String n : neighbors) {
//                    if (graph.lookupFragmentKmer(n)) {
//                        fragmentNeighbors.add(n);
//                    }
//                }
//                
//                if (fragmentNeighbors.isEmpty()) {
//                    // 0 fragment successors
//                    break;
//                }
//                else if (fragmentNeighbors.size() == 1) {
//                    // 1 fragment successor
//                    best = fragmentNeighbors.peek();
//                }
//                else {
//                    // >1 fragment successors
//                    int numKmers = kmers.size();
//                    
//                    if (numKmers >= distance) {
//                        LinkedList<String> pairedNeighbors = new LinkedList<>();
//                        
//                        partner = kmers.get(numKmers - distance);
//                        for (String n : fragmentNeighbors) {
//                            if (graph.lookupPairedKmers(n, partner)) {
//                                pairedNeighbors.add(n);
//                            }
//                        }
//                        
//                        if (pairedNeighbors.size() == 1) {
//                            best = pairedNeighbors.peek();
//                        }
//                    }
//                }
//            }
//            
//            if (best == null) {
//                break;
//            }
//            else {
//                int numKmers = kmers.size();
//                
//                if (numKmers >= distance) {
//                    if (partner == null) {
//                        partner = kmers.get(numKmers - distance);
//                    }
//                    
//                    if (transcriptKmers.contains(best)) {
//                        if (partner != null) {
//                            String pairedKmersStr = best + partner;
//                            if (!graph.lookupPairedKmers(best, partner) || usedPairs.contains(pairedKmersStr)) {
//                                break;
//                            }
//                            else {
//                                usedPairs.add(pairedKmersStr);
//                            }
//                        }
//                    }
//                }
//            }
//            
//            kmers.add(best);
//            transcriptKmers.add(best);
//            neighbors = graph.getPredecessors(best);
//        }
//        
//        Collections.reverse(kmers);
//        
//        return assembleString(kmers);
//    }
    
    private static LinkedList<Kmer> getSuccessorsRanked(Kmer source, BloomFilterDeBruijnGraph graph, int lookahead, float covThreshold) {
        LinkedList<Kmer> results = new LinkedList<>();
        ArrayDeque<Kmer> neighbors = source.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph);
        if (neighbors.isEmpty()) {
            return results;
        }
        else if (neighbors.size() == 1) {
            results.add(neighbors.peek());
            return results;
        }
        
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : neighbors) {
            if (n.count >= covThreshold) {
                float c = getMaxMedianCoverageRight(graph, n, lookahead);

                if (results.isEmpty()) {
                    results.add(n);
                    values.add(c);
                }
                else {
                    resultsItr = results.listIterator();
                    valuesItr = values.listIterator();

                    float val;
                    while (valuesItr.hasNext()) {
                        resultsItr.next();
                        val = valuesItr.next();

                        if (c > val) {
                            if (valuesItr.hasPrevious()) {
                                resultsItr.previous();
                                valuesItr.previous();
                            }

                            resultsItr.add(n);
                            valuesItr.add(c);
                            break;
                        }
                    }
                }
            }
        }
        
        return results;
    }
    
    private static LinkedList<Kmer> getPredecessorsRanked(Kmer source, BloomFilterDeBruijnGraph graph, int lookahead, float covThreshold) {
        LinkedList<Kmer> results = new LinkedList<>();
        ArrayDeque<Kmer> neighbors = source.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph);
        if (neighbors.isEmpty()) {
            return results;
        }
        else if (neighbors.size() == 1) {
            results.add(neighbors.peek());
            return results;
        }
        
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : neighbors) {
            if (n.count >= covThreshold) {
                float c = getMaxMedianCoverageLeft(graph, n, lookahead);

                if (results.isEmpty()) {
                    results.add(n);
                    values.add(c);
                }
                else {
                    resultsItr = results.listIterator();
                    valuesItr = values.listIterator();

                    float val;
                    while (valuesItr.hasNext()) {
                        resultsItr.next();
                        val = valuesItr.next();

                        if (c > val) {
                            if (valuesItr.hasPrevious()) {
                                resultsItr.previous();
                                valuesItr.previous();
                            }

                            resultsItr.add(n);
                            valuesItr.add(c);
                            break;
                        }
                    }
                }
            }
        }
        
        return results;
    }
    
    private static LinkedList<Kmer> getSuccessorsRanked(Kmer source, BloomFilterDeBruijnGraph graph, int lookahead) {
        LinkedList<Kmer> results = new LinkedList<>();
        ArrayDeque<Kmer> neighbors = source.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph);
        if (neighbors.isEmpty()) {
            return results;
        }
        else if (neighbors.size() == 1) {
            results.add(neighbors.peek());
            return results;
        }
        
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : neighbors) {
            float c = getMaxMedianCoverageRight(graph, n, lookahead);
            
            if (results.isEmpty()) {
                results.add(n);
                values.add(c);
            }
            else {
                resultsItr = results.listIterator();
                valuesItr = values.listIterator();
                
                float val;
                while (valuesItr.hasNext()) {
                    resultsItr.next();
                    val = valuesItr.next();
                    
                    if (c > val) {
                        if (valuesItr.hasPrevious()) {
                            resultsItr.previous();
                            valuesItr.previous();
                        }
                        
                        resultsItr.add(n);
                        valuesItr.add(c);
                        break;
                    }
                }
            }
        }
        
        return results;
    }
    
    private static LinkedList<Kmer> getPredecessorsRanked(Kmer source, BloomFilterDeBruijnGraph graph, int lookahead) {
        LinkedList<Kmer> results = new LinkedList<>();
        ArrayDeque<Kmer> neighbors = source.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph);
        if (neighbors.isEmpty()) {
            return results;
        }
        else if (neighbors.size() == 1) {
            results.add(neighbors.peek());
            return results;
        }
        
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : neighbors) {
            float c = getMaxMedianCoverageLeft(graph, n, lookahead);
            
            if (results.isEmpty()) {
                results.add(n);
                values.add(c);
            }
            else {
                resultsItr = results.listIterator();
                valuesItr = values.listIterator();
                
                float val;
                while (valuesItr.hasNext()) {
                    resultsItr.next();
                    val = valuesItr.next();
                    
                    if (c > val) {
                        if (valuesItr.hasPrevious()) {
                            resultsItr.previous();
                            valuesItr.previous();
                        }
                        
                        resultsItr.add(n);
                        valuesItr.add(c);
                        break;
                    }
                }
            }
        }
        
        return results;
    }
    
    private static int maxRightPartnerSearchDepth2(ArrayList<Kmer> kmers, 
                                                BloomFilterDeBruijnGraph graph,
                                                int pairedKmerDistance,
                                                BloomFilter assembledKmersBloomFilter,
                                                int minNumPairs) {
//        int k = graph.getK();
//        int numHash = graph.getMaxNumHash();
        
        final int numKmers = kmers.size();
        
//        float minCov = getMinimumKmerCoverage(kmers, Math.max(0, numKmers-pairedKmerDistance), numKmers-1);
//        int k = graph.getK();
//        int minAnchorDistanceFromEdge = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(minCov), pairedKmerDistance-k);
        
//        float[] covs = new float[Math.min(pairedKmerDistance, numKmers)];
//        for (int i=0; i<covs.length; ++i) {
//            covs[i] = kmers.get(numKmers-i-1).count;
//        }
//        
//        int k = graph.getK();
//        int minAnchorDistanceFromEdge = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(getMinium(covs)), pairedKmerDistance-k);
        
        int end = Math.max(0, numKmers-pairedKmerDistance);
        int start;

        for (start=numKmers-1; start>=end; --start) {
            if (!assembledKmersBloomFilter.lookup(kmers.get(start).getHash())) {
                break;
            }
        }
        
//        start = Math.min(start, getMinimumKmerCoverageIndexL2R(kmers, end, numKmers-1));
        
        Kmer kmer;
        Kmer p;
        end += minNumPairs;
        for (int i=start-minNumPairs; i>=end; --i) {
            kmer = kmers.get(i);
            
            if (graph.lookupLeftKmer(kmer.getHash())) {
                boolean hasComplexKmer = !graph.isLowComplexity(kmer);
                int numPaired = 1;
                for (int j=1; j<minNumPairs; ++j) {
                    p = kmers.get(i-j);
                    
                    if (graph.lookupLeftKmer(p.getHash())) {
                        if (!hasComplexKmer) {
                            hasComplexKmer = !graph.isLowComplexity(p);
                        }
                        
                        if (++numPaired >= minNumPairs && hasComplexKmer) {
                            return pairedKmerDistance - (numKmers - 1 - i);
                        }
                    }
                    else {
                        break;
                    }
                }                
            }
        }
        
//        for (int i=Math.min(numKmers-minNumPairs, start); i>=end; --i) {
//            kmer = kmers.get(i);
//            
//            if (graph.lookupLeftKmer(kmer.hashVals) &&
//                    !isLowComplexity2(kmer.bytes)) {
//                boolean found = true;
//                for (int j=1; j<minNumPairs; ++j) {
//                    p = kmers.get(i+j);
//                    if (!graph.lookupLeftKmer(p.hashVals) ||
//                            isLowComplexity2(p.bytes)) {
//                        found = false;
//                        break;
//                    }
//                }
//                if (found) {
//                    return pairedKmerDistance - (numKmers - 1 - i);
//                }
//            }
//        }
        
        return 0;
    }
    
    private static int maxLeftPartnerSearchDepth2(ArrayList<Kmer> kmers, 
                                                BloomFilterDeBruijnGraph graph,
                                                int pairedKmerDistance,
                                                BloomFilter assembledKmersBloomFilter,
                                                int minNumPairs) {
        
        // NOTE: kmers have been reversed; end of list = leftmost kmer
        
//        int k = graph.getK();
//        int numHash = graph.getMaxNumHash();
        final int numKmers = kmers.size();
        
//        float minCov = getMinimumKmerCoverage(kmers, Math.max(0, numKmers-pairedKmerDistance), numKmers-1);
//        int k = graph.getK();
//        int minAnchorDistanceFromEdge = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(minCov), pairedKmerDistance-k);
        
//        float[] covs = new float[Math.min(pairedKmerDistance, numKmers)];
//        for (int i=0; i<covs.length; ++i) {
//            covs[i] = kmers.get(numKmers-i-1).count;
//        }
//        
//        int k = graph.getK();
//        int minAnchorDistanceFromEdge = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(getMinium(covs)), pairedKmerDistance-k);
        
        int end = Math.max(0, numKmers-pairedKmerDistance);
        int start;

        for (start=numKmers-1; start>=end; --start) {
            if (!assembledKmersBloomFilter.lookup(kmers.get(start).getHash())) {
                break;
            }
        }
        
//        start = Math.min(start, getMinimumKmerCoverageIndexL2R(kmers, end, numKmers-1));
                
        Kmer kmer;
        Kmer p;
        end += minNumPairs;
        for (int i=start-minNumPairs; i>=end; --i) {
            kmer = kmers.get(i);
            
            if (graph.lookupRightKmer(kmer.getHash())) {
                boolean hasComplexKmer = !graph.isLowComplexity(kmer);
                
                int numPaired = 1;
                for (int j=1; j<minNumPairs; ++j) {
                    p = kmers.get(i-j);
                    
                    if (graph.lookupRightKmer(p.getHash())) {
                        if (!hasComplexKmer) {
                            hasComplexKmer = !graph.isLowComplexity(p);
                        }
                        
                        if (++numPaired >= minNumPairs && hasComplexKmer) {
                            return pairedKmerDistance - (numKmers - 1 - i);
                        }
                    }
                    else {
                        break;
                    }
                }                
            }
        }
        
//        for (int i=Math.min(numKmers-minNumPairs, start); i>=end; --i) {
//            kmer = kmers.get(i);
//            
//            if (graph.lookupRightKmer(kmer.hashVals) &&
//                    !isLowComplexity2(kmer.bytes)) {
//                boolean found = true;
//                for (int j=1; j<minNumPairs; ++j) {
//                    p = kmers.get(i+j);
//                    if (!graph.lookupRightKmer(p.hashVals) ||
//                            isLowComplexity2(p.bytes)) {
//                        found = false;
//                        break;
//                    }
//                }
//                if (found) {
//                    return pairedKmerDistance - (numKmers - 1 - i);
//                }
//            }
//        }
        
        return 0;
    }
    
    private static int maxRightPartnerSearchDepth2(ArrayList<Kmer> kmers, 
                                                BloomFilterDeBruijnGraph graph,
                                                int pairedKmerDistance,
                                                int minNumPairs) {
//        int k = graph.getK();
//        int numHash = graph.getMaxNumHash();
        final int numKmers = kmers.size();
                
        int end = Math.max(0, numKmers-pairedKmerDistance);
        
        Kmer kmer;
        Kmer p;
        end += minNumPairs;
        for (int i=numKmers-1-minNumPairs; i>=end; --i) {
            kmer = kmers.get(i);
            
            if (graph.lookupLeftKmer(kmer.getHash())) {
                boolean hasComplexKmer = !graph.isLowComplexity(kmer);
                
                int numPaired = 1;
                for (int j=1; j<minNumPairs; ++j) {
                    p = kmers.get(i-j);
                    
                    if (graph.lookupLeftKmer(p.getHash())) {
                        if (!hasComplexKmer) {
                            hasComplexKmer = !graph.isLowComplexity(p);
                        }
                        
                        if (++numPaired >= minNumPairs && hasComplexKmer) {
                            return pairedKmerDistance - (numKmers - 1 - i);
                        }
                    }
                    else {
                        break;
                    }
                }                
            }
        }
        
        return 0;
    }
    
    private static int maxLeftPartnerSearchDepth2(ArrayList<Kmer> kmers, 
                                                BloomFilterDeBruijnGraph graph,
                                                int pairedKmerDistance,
                                                int minNumPairs) {
        
        // NOTE: kmers have been reversed; end of list = leftmost kmer
//        int k = graph.getK();
//        int numHash = graph.getMaxNumHash();
        final int numKmers = kmers.size();
        
        int end = Math.max(0, numKmers-pairedKmerDistance);
                
        Kmer kmer;
        Kmer p;
        end += minNumPairs;
        for (int i=numKmers-1-minNumPairs; i>=end; --i) {
            kmer = kmers.get(i);
            
            if (graph.lookupRightKmer(kmer.getHash())) {
                boolean hasComplexKmer = !graph.isLowComplexity(kmer);
                
                int numPaired = 1;
                for (int j=1; j<minNumPairs; ++j) {
                    p = kmers.get(i-j);
                    
                    if (graph.lookupRightKmer(p.getHash())) {
                        if (!hasComplexKmer) {
                            hasComplexKmer = !graph.isLowComplexity(p);
                        }
                        
                        if (++numPaired >= minNumPairs && hasComplexKmer) {
                            return pairedKmerDistance - (numKmers - 1 - i);
                        }
                    }
                    else {
                        break;
                    }
                }                
            }
        }
        
        return 0;
    }
    
    
//    private static int maxRightPartnerSearchDepth(ArrayList<Kmer2> fragmentKmers, BloomFilterDeBruijnGraph graph, int pairedKmerDistance) {
//        
//        int numHash = graph.getMaxNumHash();
//        final int numKmers = fragmentKmers.size();
//        
//        float[] covs = new float[Math.min(pairedKmerDistance, numKmers)];
//        for (int i=0; i<covs.length; ++i) {
//            covs[i] = fragmentKmers.get(numKmers-i-1).count;
//        }
//        
//        int k = graph.getK();
//        int anchorLength = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(getMinium(covs)), pairedKmerDistance-k);
//        if (anchorLength == 0) {
//            ++anchorLength;
//        }
//        
//        
//        int end = Math.max(0, numKmers-pairedKmerDistance);
//        
//        Kmer2 kmer;
//        for (int i=Math.min(numKmers-anchorLength, numKmers-1); i>end; --i) {
//            kmer = fragmentKmers.get(i);
//            if (graph.lookupLeftKmer(kmer.getHashValues(k, numHash)) && i>0 && !graph.isLowComplexity(kmer)) {
//                kmer = fragmentKmers.get(i-1);
//                if (graph.lookupLeftKmer(kmer.getHashValues(k, numHash)) && !graph.isLowComplexity(kmer)) {
//                    return pairedKmerDistance - (numKmers - i);
//                }
//                else {
//                    --i;
//                }
//            }
//        }
//        
//        return 0;
//    }
//    
//    private static int maxLeftPartnerSearchDepth(ArrayList<Kmer2> fragmentKmers, BloomFilterDeBruijnGraph graph, int pairedKmerDistance, BloomFilter assembledKmersBloomFilter) {
//        
//        int numHash = graph.getMaxNumHash();
//        final int numKmers = fragmentKmers.size();
//        
//        float[] covs = new float[Math.min(pairedKmerDistance, numKmers)];
//        for (int i=0; i<covs.length; ++i) {
//            covs[i] = fragmentKmers.get(numKmers-i-1).count;
//        }
//        
//        int k = graph.getK();
//        int anchorLength = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(getMinium(covs)), pairedKmerDistance-k);
//        if (anchorLength == 0) {
//            ++anchorLength;
//        }
//        
//        int end = Math.max(0, numKmers-pairedKmerDistance);
//        int start;
//        
//        for (start=numKmers-1; start>end; --start) {
//            if (!assembledKmersBloomFilter.lookup(fragmentKmers.get(start).getHashValues(k, numHash))) {    
//                break;
//            }
//        }
//        
//        start = Math.min(numKmers-anchorLength, start);
//        
//        Kmer2 kmer;
//        for (int i=start; i>end; --i) {
//            kmer = fragmentKmers.get(i);
//            if (graph.lookupRightKmer(kmer.getHashValues(k, numHash)) && i>0 && !graph.isLowComplexity(kmer)) {
//                kmer = fragmentKmers.get(i-1);
//                if (graph.lookupRightKmer(kmer.getHashValues(k, numHash)) && !graph.isLowComplexity(kmer)) {
//                    return pairedKmerDistance - (numKmers - i);
//                }
//                else {
//                    --i;
//                }
//            }
//        }
//        
//        return 0;
//    }
    
    private static boolean hasPairedRightKmers(Kmer source,
                                            ArrayList<Kmer> kmers,
                                            int partnerFromIndex,
                                            int partnerToIndex,
                                            BloomFilterDeBruijnGraph graph) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        Kmer kmer;
        Kmer partner;
        
        partner = kmers.get(partnerFromIndex);
        
        if (graph.lookupKmerPair(partner, source) &&
                (!graph.isLowComplexity(partner) || !graph.isLowComplexity(source))) {
            
            ArrayDeque<Kmer> frontier = source.getSuccessors(k, numHash, graph);
            ArrayDeque<Kmer> newFrontier = new ArrayDeque();
            ArrayDeque<Kmer> tmp;

            Iterator<Kmer> itr;

            for (int i=partnerFromIndex+1; i<partnerToIndex; ++i) {
                partner = kmers.get(i);

                if (!graph.lookupLeftKmer(partner.getHash())) {
                    return false;
                }

                boolean partnerIsLowComplexity = graph.isLowComplexity(partner);
                
                itr = frontier.iterator();
                while (itr.hasNext()) {
                    kmer = itr.next();
                    if (graph.lookupRightKmer(kmer.getHash()) && 
                            graph.lookupKmerPairing(partner, kmer) &&
                            (!partnerIsLowComplexity || !graph.isLowComplexity(kmer))) {
                        newFrontier.addAll(kmer.getSuccessors(k, numHash, graph));
                    }

                    itr.remove();
                }

                if (newFrontier.isEmpty()) {
                    return false;
                }

                tmp = frontier; // an empty list
                frontier = newFrontier; // neighbors
                newFrontier = tmp; // reuse empty list
            }
        }
        else {
            return false;
        }
        
        return true;
    }
   
    private static boolean hasPairedLeftKmers(Kmer source,
                                            ArrayList<Kmer> kmers,
                                            int partnerFromIndex,
                                            int partnerToIndex,
                                            BloomFilterDeBruijnGraph graph) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        Kmer kmer;
        Kmer partner;
        
        partner = kmers.get(partnerFromIndex);
        
        if (graph.lookupKmerPair(source, partner) &&
                (!graph.isLowComplexity(source) || !graph.isLowComplexity(partner))) {
            
            ArrayDeque<Kmer> frontier = source.getPredecessors(k, numHash, graph);
            ArrayDeque<Kmer> newFrontier = new ArrayDeque<>();
            ArrayDeque<Kmer> tmp;

            Iterator<Kmer> itr;

            for (int i=partnerFromIndex+1; i<partnerToIndex; ++i) {
                partner = kmers.get(i);

                if (!graph.lookupRightKmer(partner.getHash())) {
                    return false;
                }
                
                boolean partnerIsLowComplexity = graph.isLowComplexity(partner);

                itr = frontier.iterator();
                while (itr.hasNext()) {
                    kmer = itr.next();
                    if (graph.lookupLeftKmer(kmer.getHash()) && 
                            graph.lookupKmerPairing(kmer, partner) &&
                            (!partnerIsLowComplexity || !graph.isLowComplexity(kmer))) {
                        newFrontier.addAll(kmer.getPredecessors(k, numHash, graph));
                    }

                    itr.remove();
                }

                if (newFrontier.isEmpty()) {
                    return false;
                }

                tmp = frontier; // an empty list
                frontier = newFrontier; // neighbors
                newFrontier = tmp; // reuse empty list
            }
        }
        else {
            return false;
        }
        
        return true;
    }
    
    private static boolean areLeftKmers(ArrayList<Kmer> kmers, 
                                        int fromIndex,
                                        int toIndex,
                                        BloomFilterDeBruijnGraph graph) {        
        for (int i=fromIndex; i<toIndex; ++i) {
            if (!graph.lookupLeftKmer(kmers.get(i).getHash())) {
                return false;
            }
        }
        
        return true;
    }
    
    private static boolean areRightKmers(ArrayList<Kmer> kmers,
                                        int fromIndex,
                                        int toIndex,
                                        BloomFilterDeBruijnGraph graph) {
        for (int i=fromIndex; i<toIndex; ++i) {
            if (!graph.lookupRightKmer(kmers.get(i).getHash())) {
                return false;
            }
        }
        
        return true;
    }
    
    private static boolean extendRightWithPairedKmersOnly(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            HashSet<Kmer> usedKmers,
                                            BloomFilter assembledKmersBloomFilter) {
        // each extension must be a paired kmer        
        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        int distance = graph.getFragPairedKmerDistance();
        int numKmers = kmers.size();
        
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = cursor.getSuccessors(k, numHash, graph);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        int partnerIndex = numKmers-distance;
        
        while (!neighbors.isEmpty()) {            
            if (numKmers < distance) {
                // too short to have partners
                return true;
            }
            
            cursor = null;
                        
            if (neighbors.isEmpty()) {
                return false;
            }
            else {
                Kmer partner = kmers.get(partnerIndex);
                            
                if (!graph.lookupLeftKmer(partner.getHash())) {
                    // not a supporting partner
                    return true;
                }
                
                for (Kmer n : neighbors) {
                    if (graph.lookupRightKmer(n.getHash()) && graph.lookupKmerPairing(partner, n)) {
                        if (cursor == null) {
                            cursor = n;
                        }
                        else {
                            // too many candidates
                            return true;
                        }
                    }
                }
                 
                if (cursor == null) {
                    // no good candidates
                    return true;
                }
                
                if (usedKmers.contains(cursor)) {
                    // check whether this kmer pair has been used in this sequence already
                    int i = kmers.size();
                    while (--i >= distance) {
                        if (cursor.equals(kmers.get(i)) &&
                                partner.equals(kmers.get(i-distance))) {
                            return false;
                        }
                    }
                }
                
                if (assembledKmersBloomFilter.lookup(cursor.getHash())) {
                    for (int i=partnerIndex; i<numKmers; ++i) {
                        if (assembledKmersBloomFilter.lookup(kmers.get(i).getHash())) {
                            if (i == numKmers-1) {
                                return false;
                            }
                        }
                        else {
                            break;
                        }
                    }
                }
                
                kmers.add(cursor);
                usedKmers.add(cursor);

                ++partnerIndex;
                ++numKmers;

                neighbors = cursor.getSuccessors(k, numHash, graph);
            }
        }
        
        return false;
    }
    
    private static ArrayDeque<Kmer> extendRightWithPairedKmersOnly(Kmer kmer,
                                                        ArrayList<Kmer> kmers, 
                                                        BloomFilterDeBruijnGraph graph,
                                                        int maxLength) {
        // each extension must be a paired kmer
        
        final int distance = graph.getFragPairedKmerDistance();
        final int numKmers = kmers.size();
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        
        if (numKmers < distance) {
            return extension;
        }
        
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();

        ArrayDeque<Kmer> neighbors = kmer.getSuccessors(k, numHash, graph);
        if (!neighbors.isEmpty()) {
            int partnerIndex = numKmers-distance+1;

            while (!neighbors.isEmpty()) {            
                if (partnerIndex >= numKmers || extension.size() >= maxLength) {
                    break;
                }

                if (neighbors.isEmpty()) {
                    break;
                }
                else {          
                    if (!graph.lookupLeftKmer(kmers.get(partnerIndex).getHash())) {
                        // not a supporting partner
                        return extension;
                    }

                    Iterator<Kmer> itr = neighbors.iterator();
                    while (itr.hasNext()) {
                        if (!graph.lookupRightKmer(itr.next().getHash())) {
                            itr.remove();
                        }
                    }

                    if (neighbors.isEmpty()) {
                        break;
                    }
                    else if (neighbors.size() == 1) {
                        extension.add(neighbors.pop());
                        ++partnerIndex;
                    }
                    else {
                        itr = neighbors.iterator();

                        ArrayList<Kmer> f = new ArrayList<>(distance);
                        for (int i=partnerIndex; i<numKmers; ++i) {
                            f.add(kmers.get(i));
                        }
                        f.addAll(extension);

                        Kmer best = itr.next();
                        ArrayDeque<Kmer> bestExtension = extendRightWithPairedKmersOnly(best,
                                                                            f, 
                                                                            graph,
                                                                            maxLength - 1 - extension.size());

                        while (itr.hasNext()) {
                            Kmer n = itr.next();
                            
                            ArrayDeque<Kmer> e = extendRightWithPairedKmersOnly(n,
                                                                            f, 
                                                                            graph,
                                                                            maxLength - 1 - extension.size());
                            
                            if (e.size() > bestExtension.size()) {
                                best = n;
                                bestExtension = e;
                            }
                        }
                        
                        extension.add(best);
                        extension.addAll(bestExtension);
                        partnerIndex += bestExtension.size();
                    }

                    neighbors = extension.getLast().getSuccessors(k, numHash, graph);
                }
            }
        }
        
        return extension;
    }
    
    private static ArrayDeque<Kmer> extendLeftWithPairedKmersOnly(Kmer kmer,
                                                        ArrayList<Kmer> kmers, 
                                                        BloomFilterDeBruijnGraph graph,
                                                        int maxLength) {
        // each extension must be a paired kmer
        
        final int distance = graph.getFragPairedKmerDistance();
        final int numKmers = kmers.size();
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        
        if (numKmers < distance) {
            return extension;
        }
        
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();

        ArrayDeque<Kmer> neighbors = kmer.getPredecessors(k, numHash, graph);
        if (!neighbors.isEmpty()) {
            int partnerIndex = numKmers-distance+1;

            while (!neighbors.isEmpty()) {            
                if (partnerIndex >= numKmers || extension.size() >= maxLength) {
                    break;
                }

                if (neighbors.isEmpty()) {
                    break;
                }
                else {          
                    if (!graph.lookupRightKmer(kmers.get(partnerIndex).getHash())) {
                        // not a supporting partner
                        return extension;
                    }

                    Iterator<Kmer> itr = neighbors.iterator();
                    while (itr.hasNext()) {
                        if (!graph.lookupLeftKmer(itr.next().getHash())) {
                            itr.remove();
                        }
                    }

                    if (neighbors.isEmpty()) {
                        break;
                    }
                    else if (neighbors.size() == 1) {
                        extension.add(neighbors.pop());
                        ++partnerIndex;
                    }
                    else {
                        itr = neighbors.iterator();

                        ArrayList<Kmer> f = new ArrayList<>(distance);
                        for (int i=partnerIndex; i<numKmers; ++i) {
                            f.add(kmers.get(i));
                        }
                        f.addAll(extension);

                        Kmer best = itr.next();
                        ArrayDeque<Kmer> bestExtension = extendLeftWithPairedKmersOnly(best,
                                                                            f, 
                                                                            graph,
                                                                            maxLength - 1 - extension.size());

                        while (itr.hasNext()) {
                            Kmer n = itr.next();
                            
                            ArrayDeque<Kmer> e = extendLeftWithPairedKmersOnly(n,
                                                                            f, 
                                                                            graph,
                                                                            maxLength - 1 - extension.size());
                            
                            if (e.size() > bestExtension.size()) {
                                best = n;
                                bestExtension = e;
                            }
                        }
                        
                        extension.add(best);
                        extension.addAll(bestExtension);
                        partnerIndex += bestExtension.size();
                    }

                    neighbors = extension.getLast().getPredecessors(k, numHash, graph);
                }
            }
        }
        
        return extension;
    }
    
    private static boolean extendRightWithPairedKmersBFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            HashSet<Kmer> usedKmers,
                                            float minKmerCov) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        int distance = graph.getFragPairedKmerDistance();
        int numKmers = kmers.size();
                            
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = cursor.getSuccessors(k, numHash, graph);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        int partnerIndex = numKmers-distance;
        
        while (!neighbors.isEmpty()) {
            ArrayDeque<Kmer> simpleExtension = extendRight(cursor,
                                            graph, 
                                            maxTipLength, 
                                            usedKmers, 
                                            maxIndelSize, 
                                            percentIdentity,
                                            minKmerCov);
                        
            if (!simpleExtension.isEmpty()) {
                kmers.addAll(simpleExtension);
                
                numKmers = kmers.size();

                if (numKmers < distance) {
                    // too short to have partners
                    return true;
                }

                partnerIndex += simpleExtension.size();

                neighbors = simpleExtension.getLast().getSuccessors(k, numHash, graph);
                if (neighbors.isEmpty()) {
                    return false;
                }
            }
            
            if (numKmers < distance) {
                // too short to have partners
                return true;
            }
            
            cursor = null;
                        
            if (neighbors.isEmpty()) {
                return false;
            }
            else if (neighbors.size() == 1) {
                cursor = neighbors.pop();
            }
            else {
                if (!areLeftKmers(kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                    // not enough supporting partners
                    return true;
                }
                                
                for (Kmer n : neighbors) {
                    if (hasPairedRightKmers(n, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                        if (cursor == null) {
                            cursor = n;
                        }
                        else {
                            return true;
                        }
                    }
                }
                              
                if (cursor == null) {
                    // no good candidates
                    return true;
                }
            }
            
            Kmer partner = kmers.get(partnerIndex);
            
            if (usedKmers.contains(cursor)) {
                int i = kmers.size();
                while (--i >= distance) {
                    if (cursor.equals(kmers.get(i)) &&
                            partner.equals(kmers.get(i-distance))) {
                        return false;
                    }
                }
            }
            
            kmers.add(cursor);
            usedKmers.add(cursor);
            
            ++partnerIndex;
            ++numKmers;
            
            neighbors = cursor.getSuccessors(k, numHash, graph);
        }
        
        return false;
    }
    
    private static boolean extendLeftWithPairedKmersOnly(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            HashSet<Kmer> usedKmers,
                                            BloomFilter assembledKmersBloomFilter) {
        // each extension must be a paired kmer        
        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        int distance = graph.getFragPairedKmerDistance();
        int numKmers = kmers.size();
        
        // Note that `kmers` are in reverse order already
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = cursor.getPredecessors(k, numHash, graph);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        int partnerIndex = numKmers-distance;
        
        while (!neighbors.isEmpty()) {            
            if (numKmers < distance) {
                // too short to have partners
                return true;
            }
            
            cursor = null;
                        
            if (neighbors.isEmpty()) {
                return false;
            }
            else {
                Kmer partner = kmers.get(partnerIndex);
                            
                if (!graph.lookupRightKmer(partner.getHash())) {
                    // not a supporting partner
                    return true;
                }
                
                for (Kmer n : neighbors) {
                    if (graph.lookupLeftKmer(n.getHash()) && graph.lookupKmerPairing(n, partner)) {
                        if (cursor == null) {
                            cursor = n;
                        }
                        else {
                            return true;
                        }
                    }
                }
                 
                if (cursor == null) {
                    // no good candidates
                    return true;
                }
                
                if (usedKmers.contains(cursor)) {
                    // check whether this kmer pair has been used in this sequence already (ie. a loop)
                    int i = kmers.size();
                    while (--i >= distance) {
                        if (cursor.equals(kmers.get(i)) &&
                                partner.equals(kmers.get(i-distance))) {
                            return false;
                        }
                    }
                }
                
                if (assembledKmersBloomFilter.lookup(cursor.getHash())) {
                    for (int i=partnerIndex; i<numKmers; ++i) {
                        if (assembledKmersBloomFilter.lookup(kmers.get(i).getHash())) {
                            if (i == numKmers-1) {
                                return false;
                            }
                        }
                        else {
                            break;
                        }
                    }
                }
                
                kmers.add(cursor);
                usedKmers.add(cursor);

                ++partnerIndex;
                ++numKmers;

                neighbors = cursor.getPredecessors(k, numHash, graph);
            }
        }
        
        return false;
    }
    
    private static boolean extendLeftWithPairedKmersBFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            HashSet<Kmer> usedKmers,
                                            float minKmerCov) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        int distance = graph.getFragPairedKmerDistance();
        int numKmers = kmers.size();
        
        // Note that `kmers` are in reverse order already
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = cursor.getPredecessors(k, numHash, graph);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        int partnerIndex = numKmers-distance;
        
        while (!neighbors.isEmpty()) {
            ArrayDeque<Kmer> simpleExtension = extendLeft(cursor,
                                            graph, 
                                            maxTipLength, 
                                            usedKmers, 
                                            maxIndelSize, 
                                            percentIdentity,
                                            minKmerCov);
 
            if (!simpleExtension.isEmpty()) {
                kmers.addAll(simpleExtension);
                
                numKmers = kmers.size();

                if (numKmers < distance) {
                    // too short to have partners
                    return true;
                }

                partnerIndex += simpleExtension.size();

                neighbors = simpleExtension.getLast().getPredecessors(k, numHash, graph);
                if (neighbors.isEmpty()) {
                    return false;
                }
            }
            
            if (numKmers < distance) {
                // too short to have partners
                return true;
            }
            
            cursor = null;
                        
            if (neighbors.isEmpty()) {
                return false;
            }
            else if (neighbors.size() == 1) {
                cursor = neighbors.pop();
            }
            else {
                if (!areRightKmers(kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                    // not enough supporting partners
                    return true;
                }
                
                for (Kmer n : neighbors) {
                    if (hasPairedLeftKmers(n, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                        if (cursor == null) {
                            cursor = n;
                        }
                        else {
                            return true;
                        }
                    }
                }
                 
                if (cursor == null) {
                    // no good candidates
                    return true;
                }
            }
            
            Kmer partner = kmers.get(partnerIndex);
            
            if (usedKmers.contains(cursor)) {
                int i = kmers.size();
                while (--i >= distance) {
                    if (cursor.equals(kmers.get(i)) &&
                            partner.equals(kmers.get(i-distance))) {
                        return false;
                    }
                }
            }
                        
            kmers.add(cursor);
            usedKmers.add(cursor);
                        
            ++partnerIndex;
            ++numKmers;
            
            neighbors = cursor.getPredecessors(k, numHash, graph);
        }
        
        return false;
    }
    
    private static boolean extendRightWithPairedKmersBFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            BloomFilter assembledKmersBloomFilter,
                                            HashSet<Kmer> usedKmers,
                                            float maxCovGradient,
                                            float minKmerCov) {
        
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int distance = graph.getFragPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int numKmers = kmers.size();
                            
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = cursor.getSuccessors(k, numHash, graph);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        int partnerIndex = numKmers-distance;
//        final int minFragmentCoverageThreshold = (int) Math.floor(1.0/maxCovGradient);
        
        while (!neighbors.isEmpty()) {
            ArrayDeque<Kmer> simpleExtension = extendRight(cursor,
                                            graph, 
                                            maxTipLength, 
                                            usedKmers, 
                                            maxIndelSize, 
                                            percentIdentity,
                                            minKmerCov);
                        
            if (!simpleExtension.isEmpty()) {
//                Iterator<Kmer> itr = simpleExtension.descendingIterator();
//                while(itr.hasNext() && assembledKmersBloomFilter.lookup(itr.next().getHash())) {
//                    itr.remove();
//                }
////                
//                if (!simpleExtension.isEmpty()) {
                    kmers.addAll(simpleExtension);
                    usedKmers.addAll(simpleExtension);
                    
                    numKmers = kmers.size();

                    if (numKmers < distance) {
                        // too short to have partners
                        return true;
                    }

                    partnerIndex += simpleExtension.size();

                    // NOTE: kmer at `partnerIndex` will be paired with `cursor`
                    
                    neighbors = simpleExtension.getLast().getSuccessors(k, numHash, graph);
                    if (neighbors.isEmpty()) {
                        return false;
                    }
//                }
            }
            
            if (numKmers < distance) {
                // too short to have partners
                return true;
            }
            
            cursor = null;
                        
            
            if (neighbors.isEmpty()) {
                return false;
            }
            else if (neighbors.size() == 1) {
                cursor = neighbors.pop();
            }
            else {
                if (!areLeftKmers(kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                    // not enough supporting partners
                    return true;
                }
                
                Iterator<Kmer> itr = neighbors.iterator();
                while (itr.hasNext()) {
                    Kmer kmer = itr.next();
                    
                    if (!hasPairedRightKmers(kmer, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                        itr.remove();
                    }
                }
                
                if (neighbors.isEmpty()) {
                    // no good candidates
                    return true;
                }
                else if (neighbors.size() == 1) {
                    cursor = neighbors.pop();
                }
                else {
                    // two or more neighbors are supported by paired kmers
                    ArrayDeque<Kmer> e = extendRightPE(kmers, graph, maxTipLength, minKmerCov);
                    if (e == null || e.isEmpty()) {
                        return false; // ambiguous branches supported by paired kmers
                    }
                    else {
                        cursor = e.getFirst();
                    }
                    
//                    float medEdgeCoverage = getMedianKmerCoverage(kmers, numKmers-1-lookahead, numKmers-1);
//                    float minFragmentCoverage = getMinimumKmerCoverage(kmers, numKmers-distance, numKmers-1);
////                    if (minFragmentCoverage >= minFragmentCoverageThreshold) {
//                        ArrayDeque<Kmer> neighborsBackUp = new ArrayDeque<>(neighbors);
//                        
//                        // only remove neighbors based on coverage when fragment coverage is not too low
//                        float minCovThreshold = minFragmentCoverage * maxCovGradient;            
//                        itr = neighbors.iterator();
//                        while (itr.hasNext()) {
//                            Kmer kmer = itr.next();
//
//                            if (kmer.count <= minCovThreshold) {
//                                itr.remove();
//                            }
//                        }
//                        
//                        if (neighbors.size() > 1 && medEdgeCoverage > minFragmentCoverage) {
//                            minCovThreshold = medEdgeCoverage * maxCovGradient;
//                            
//                            itr = neighbors.iterator();
//                            while (itr.hasNext()) {
//                                Kmer kmer = itr.next();
//
//                                if (kmer.count <= minCovThreshold) {
//                                    itr.remove();
//                                }
//                            }
//                        }
//
//                        if (neighbors.size() == 1) {
//                            cursor = neighbors.pop();
//                        }
//                        else {
//                            if (neighbors.isEmpty()) {
//                                neighbors = neighborsBackUp;
//                            }
//                            
//                            Kmer best = null;
//                            float bestCov = 0;
//                            float secondCov = 0;
//                            ArrayDeque<Kmer> bestGE = null;
//                            ArrayDeque<Kmer> secondGE = null;
//                            
//                            for (Kmer kmer : neighbors) {
//                                
//                                ArrayDeque<Kmer> e = extendRightWithPairedKmersOnly(kmer, kmers, graph, lookahead);
//                                e.addFirst(kmer);
//                                
//                                ArrayDeque<Kmer> ge = greedyExtendRight(graph, e.getLast(), lookahead, k-lookahead);
//                                Iterator<Kmer> ditr = e.descendingIterator();
//                                while (ditr.hasNext()) {
//                                    ge.addFirst(ditr.next());
//                                }
//                                
//                                float c = getMedianKmerCoverage(e);
//                                
//                                if (best == null) {
//                                    best = kmer;
//                                    bestCov = c;
//                                    bestGE = ge;
//                                } 
//                                else {
//                                    if (bestCov < c) {
//                                        secondCov = bestCov;
//                                        secondGE = bestGE;
//                                        
//                                        best = kmer;
//                                        bestCov = c;
//                                        bestGE = ge;
//                                    }
//                                    else if (secondCov < c) {
//                                        secondCov = c;
//                                        secondGE = ge;
//                                    }
//                                }
//                            }
//
//                            if (bestCov * maxCovGradient >= secondCov || 
//                                    secondCov - bestCov * maxCovGradient <= 0.1f * secondCov ||
//                                    getPercentIdentity(graph.assemble(bestGE), graph.assemble(secondGE)) >= percentIdentity) {
//                                cursor = best;
//                            }
//                            else {
//                                ArrayDeque<Kmer> e = extendRightPE(kmers, graph, maxTipLength, maxIndelSize, percentIdentity);
//                                if (e == null || e.isEmpty()) {
//                                    return false; // ambiguous branches supported by paired kmers
//                                }
//                                else {
//                                    cursor = e.getLast();
//                                }
//                            }
//                        }
////                    }
////                    else {
////                        return false; // ambiguous branches supported by paired kmers
////                    }
                }
            }
            
            if (usedKmers.contains(cursor)){
                // check whether this kmer pair has been used in this sequence already
                Kmer partner = kmers.get(partnerIndex);
                for (int i = kmers.size()-1; i >= distance; --i) {
                    if (cursor.equals(kmers.get(i)) &&
                            partner.equals(kmers.get(i-distance))) {
                        // kmer pair has been used previously
                        return false;
                    }
                }
            }
            
            if (assembledKmersBloomFilter.lookup(cursor.getHash())) {
                boolean assembled = greedyExtendRight(graph, cursor, lookahead, lookahead, assembledKmersBloomFilter) != null;
                
                if (assembled) {
//                    int numNotAssembled = 0;
                    for (int i=partnerIndex; i<numKmers; ++i) {
                        if (!assembledKmersBloomFilter.lookup(kmers.get(i).getHash())) {
//                            if (distanceInversePI < ++numNotAssembled) {
                                assembled = false;
                                break;
//                            }
                        }
                    }
                }
                
                if (assembled) {
                    // region already assembled, do not extend further
                    return false;
                }
            }
            
            kmers.add(cursor);
            
            usedKmers.add(cursor);
            
            ++partnerIndex;
            ++numKmers;
            
            neighbors = cursor.getSuccessors(k, numHash, graph);
        }
        
        return false;
    }
        
    private static boolean extendLeftWithPairedKmersBFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            BloomFilter assembledKmersBloomFilter,
                                            HashSet<Kmer> usedKmers,
                                            float maxCovGradient,
                                            float minKmerCov) {
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        
        final int distance = graph.getFragPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int numKmers = kmers.size();
        
        // Note that `kmers` are in reverse order already
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = cursor.getPredecessors(k, numHash, graph);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        int partnerIndex = numKmers-distance;
//        final int minFragmentCoverageThreshold = (int) Math.floor(1.0/maxCovGradient);
        
        while (!neighbors.isEmpty()) {
            ArrayDeque<Kmer> simpleExtension = extendLeft(cursor,
                                            graph, 
                                            maxTipLength, 
                                            usedKmers, 
                                            maxIndelSize, 
                                            percentIdentity,
                                            minKmerCov);
            
            if (!simpleExtension.isEmpty()) {
//                Iterator<Kmer> itr = simpleExtension.descendingIterator();
//                while(itr.hasNext() && assembledKmersBloomFilter.lookup(itr.next().getHash())) {
//                    itr.remove();
//                }
//                
//                if (!simpleExtension.isEmpty()) {
                    kmers.addAll(simpleExtension);
                    usedKmers.addAll(simpleExtension);
                    
                    numKmers = kmers.size();

                    if (numKmers < distance) {
                        // too short to have partners
                        return true;
                    }

                    partnerIndex += simpleExtension.size();

                    // NOTE: kmer at `partnerIndex` will be paired with `cursor`
                    
                    neighbors = simpleExtension.getLast().getPredecessors(k, numHash, graph);
                    if (neighbors.isEmpty()) {
                        return false;
                    }
//                }
            }
            
            if (numKmers < distance) {
                // too short to have partners
                return true;
            }
            
            cursor = null;
                                    
            if (neighbors.isEmpty()) {
                return false;
            }
            else if (neighbors.size() == 1) {
                cursor = neighbors.pop();
            }
            else {
                if (!areRightKmers(kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                    // not enough supporting partners
                    return true;
                }
                
                Iterator<Kmer> itr = neighbors.iterator();
                while (itr.hasNext()) {
                    Kmer kmer = itr.next();
                    
                    if (!hasPairedLeftKmers(kmer, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                        itr.remove();
                    }
                }
                
                if (neighbors.isEmpty()) {
                    // no good candidates
                    return true;
                }
                else if (neighbors.size() == 1) {
                    cursor = neighbors.pop();
                }
                else {
                    // two or more neighbors are supported by paired kmers
                    
                    ArrayDeque<Kmer> e = extendLeftPE(kmers, graph, maxTipLength, minKmerCov);
                    if (e == null || e.isEmpty()) {
                        return false; // ambiguous branches supported by paired kmers
                    }
                    else {
                        cursor = e.getFirst();
                    }
                    
//                    float medEdgeCoverage = getMedianKmerCoverage(kmers, numKmers-1-lookahead, numKmers-1);
//                    float minFragmentCoverage = getMinimumKmerCoverage(kmers, numKmers-distance, numKmers-1);
////                    if (minFragmentCoverage >= minFragmentCoverageThreshold) {
//                        ArrayDeque<Kmer> neighborsBackUp = new ArrayDeque<>(neighbors);
//                        
//                        // only remove neighbors based on coverage when fragment coverage is not too low
//                        float minCovThreshold = minFragmentCoverage * maxCovGradient;            
//                        itr = neighbors.iterator();
//                        while (itr.hasNext()) {
//                            Kmer kmer = itr.next();
//
//                            if (kmer.count <= minCovThreshold) {
//                                itr.remove();
//                            }
//                        }
//                        
//                        if (neighbors.size() > 1 && medEdgeCoverage > minFragmentCoverage) {
//                            minCovThreshold = medEdgeCoverage * maxCovGradient;
//                            
//                            itr = neighbors.iterator();
//                            while (itr.hasNext()) {
//                                Kmer kmer = itr.next();
//
//                                if (kmer.count <= minCovThreshold) {
//                                    itr.remove();
//                                }
//                            }
//                        }
//                        
//                        if (neighbors.size() == 1) {
//                            cursor = neighbors.pop();
//                        }
//                        else {
//                            if (neighbors.isEmpty()) {
//                                neighbors = neighborsBackUp;
//                            }
//                            
//                            Kmer best = null;
//                            float bestCov = 0;
//                            float secondCov = 0;
//                            ArrayDeque<Kmer> bestGE = null;
//                            ArrayDeque<Kmer> secondGE = null;
//                            
//                            for (Kmer kmer : neighbors) {
//                                
//                                ArrayDeque<Kmer> e = extendLeftWithPairedKmersOnly(kmer, kmers, graph, lookahead);
//                                e.addFirst(kmer);
//                                
//                                ArrayDeque<Kmer> ge = greedyExtendLeft(graph, e.getLast(), lookahead, k-lookahead);
//                                Iterator<Kmer> ditr = e.descendingIterator();
//                                while (ditr.hasNext()) {
//                                    ge.addFirst(ditr.next());
//                                }
//                                
//                                float c = getMedianKmerCoverage(e);
//                                
//                                if (best == null) {
//                                    best = kmer;
//                                    bestCov = c;
//                                    bestGE = ge;
//                                }
//                                else {
//                                    if (bestCov < c) {
//                                        secondCov = bestCov;
//                                        secondGE = bestGE;
//                                        
//                                        best = kmer;
//                                        bestCov = c;
//                                        bestGE = ge;
//                                    }
//                                    else if (secondCov < c) {
//                                        secondCov = c;
//                                        secondGE = ge;
//                                    }
//                                }
//                            }
//                            
//                            if (bestCov * maxCovGradient >= secondCov || 
//                                    secondCov - bestCov * maxCovGradient <= 0.1f * secondCov ||
//                                    getPercentIdentity(graph.assembleReverseOrder(bestGE), graph.assembleReverseOrder(secondGE)) >= percentIdentity) {
//                                cursor = best;
//                            }
//                            else {
//                                return false; // ambiguous branches supported by paired kmers
//                            }
//                        }
////                    }
////                    else {
////                        return false; // ambiguous branches supported by paired kmers
////                    }
                }
            }
            
            if (usedKmers.contains(cursor)){
                // check whether this kmer pair has been used in this sequence already (ie. a loop)
                Kmer partner = kmers.get(partnerIndex);
                for (int i = kmers.size()-1; i >= distance; --i) {
                    if (cursor.equals(kmers.get(i)) &&
                            partner.equals(kmers.get(i-distance))) {
                        return false;
                    }
                }
            }
            
            if (assembledKmersBloomFilter.lookup(cursor.getHash())) {
                boolean assembled = greedyExtendLeft(graph, cursor, lookahead, lookahead, assembledKmersBloomFilter) != null;
                
                if (assembled) {
//                    int numNotAssembled = 0;
                    for (int i=partnerIndex; i<numKmers; ++i) {
                        if (!assembledKmersBloomFilter.lookup(kmers.get(i).getHash())) {
//                            if (distanceInversePI < ++numNotAssembled) {
                                assembled = false;
                                break;
//                            }
                        }
                    }
                }
                
                if (assembled) {
                    // region already assembled, do not extend further
                    return false;
                }
            }
            
            kmers.add(cursor);
            
            usedKmers.add(cursor);
                        
            ++partnerIndex;
            ++numKmers;
            
            neighbors = cursor.getPredecessors(k, numHash, graph);
        }
        
        return false;
    }
    
    private static void countKmerPairsSE(final ArrayList<Kmer> leftKmers,
                                    final ArrayDeque<Kmer> rightKmers,
                                    final int gap,
                                    final BloomFilterDeBruijnGraph graph,
                                    final int[] result) {

        final int readPairedKmersDist = graph.getReadPairedKmerDistance();

        int supportingReadKmerPairs = 0;
        
        final int numLeftKmers = leftKmers.size();
        final int maxRightKmersIndex = Math.min(readPairedKmersDist - 1 - gap, rightKmers.size() - 1);
        
        int leftReadPartnerIndex = numLeftKmers - readPairedKmersDist + gap;
        int lastSupportedKmerIndex = -1;
        
        Iterator<Kmer> itr = rightKmers.iterator();
        for (int i=0; i<=maxRightKmersIndex; ++i) {
            Kmer rightKmer = itr.next();
            if (leftReadPartnerIndex >=0 && leftReadPartnerIndex < numLeftKmers) {
                if (graph.lookupReadKmerPair(leftKmers.get(leftReadPartnerIndex), rightKmer)) {
                    ++supportingReadKmerPairs;
                    lastSupportedKmerIndex = i;
                }
            }
            
            ++leftReadPartnerIndex;
            
            if (leftReadPartnerIndex >= numLeftKmers) {
                break;
            }
        }
        
        result[0] = supportingReadKmerPairs;
        result[1] = lastSupportedKmerIndex;
    }
    
    private static void countKmerPairsReversedSE(final ArrayDeque<Kmer> leftKmers,
                                    final ArrayList<Kmer> rightKmers,
                                    final int gap,
                                    final BloomFilterDeBruijnGraph graph,
                                    final int[] result) {

        final int readPairedKmersDist = graph.getReadPairedKmerDistance();

        int supportingReadKmerPairs = 0;
        
        final int numRightKmers = rightKmers.size();
        final int maxLeftKmersIndex = Math.min(readPairedKmersDist - 1 - gap, leftKmers.size() - 1);
        
        int rightReadPartnerIndex = numRightKmers - readPairedKmersDist + gap;
        int lastSupportedKmerIndex = -1;
        
        Iterator<Kmer> itr = leftKmers.iterator();
        for (int i=0; i<=maxLeftKmersIndex; ++i) {
            Kmer leftKmer = itr.next();
            if (rightReadPartnerIndex >=0 && rightReadPartnerIndex < numRightKmers) {
                if (graph.lookupReadKmerPair(leftKmer, rightKmers.get(rightReadPartnerIndex))) {
                    ++supportingReadKmerPairs;
                    lastSupportedKmerIndex = i;
                }
            }
            
            ++rightReadPartnerIndex;
            
            if (rightReadPartnerIndex >= numRightKmers) {
                break;
            }
        }
        
        result[0] = supportingReadKmerPairs;
        result[1] = lastSupportedKmerIndex;
    }
    
    private static void countKmerPairsPE(final ArrayList<Kmer> leftKmers,
                                    final ArrayDeque<Kmer> rightKmers,
                                    final int gap,
                                    final BloomFilterDeBruijnGraph graph,
                                    final int[] result) {

        final int readPairedKmersDist = graph.getReadPairedKmerDistance();
        final int fragPairedKmersDist = graph.getFragPairedKmerDistance();

        int supportingReadKmerPairs = 0;
        int supportingFragKmerPairs = 0;
        
        final int numLeftKmers = leftKmers.size();
        final int maxRightKmersIndex = Math.min(fragPairedKmersDist - 1 - gap, rightKmers.size() - 1);
        
        int leftReadPartnerIndex = numLeftKmers - readPairedKmersDist + gap;
        int leftFragPartnerIndex = numLeftKmers - fragPairedKmersDist + gap;
        int lastSupportedKmerIndex = -1;
        
        Iterator<Kmer> itr = rightKmers.iterator();
        for (int i=0; i<=maxRightKmersIndex; ++i) {
            Kmer rightKmer = itr.next();
            if (leftReadPartnerIndex >=0 && leftReadPartnerIndex < numLeftKmers) {
                if (graph.lookupReadKmerPair(leftKmers.get(leftReadPartnerIndex), rightKmer)) {
                    ++supportingReadKmerPairs;
                    lastSupportedKmerIndex = i;
                }
            }
            
            if (leftFragPartnerIndex >= 0 && leftFragPartnerIndex < numLeftKmers) {
                if (graph.lookupKmerPair(leftKmers.get(leftFragPartnerIndex), rightKmer)) {
                    ++supportingFragKmerPairs;
                    lastSupportedKmerIndex = i;
                }
            }
            
            ++leftReadPartnerIndex;
            ++leftFragPartnerIndex;
            
            if (leftReadPartnerIndex >= numLeftKmers && leftFragPartnerIndex >= numLeftKmers) {
                break;
            }
        }
        
        result[0] = supportingReadKmerPairs;
        result[1] = supportingFragKmerPairs;
        result[2] = lastSupportedKmerIndex;
    }
    
    private static void countKmerPairsReversedPE(final ArrayDeque<Kmer> leftKmers,
                                    final ArrayList<Kmer> rightKmers,
                                    final int gap,
                                    final BloomFilterDeBruijnGraph graph,
                                    final int[] result) {

        final int readPairedKmersDist = graph.getReadPairedKmerDistance();
        final int fragPairedKmersDist = graph.getFragPairedKmerDistance();

        int supportingReadKmerPairs = 0;
        int supportingFragKmerPairs = 0;
        
        final int numRightKmers = rightKmers.size();
        final int maxLeftKmersIndex = Math.min(fragPairedKmersDist - 1 - gap, leftKmers.size() - 1);
        
        int rightReadPartnerIndex = numRightKmers - readPairedKmersDist + gap;
        int rightFragPartnerIndex = numRightKmers - fragPairedKmersDist + gap;
        int lastSupportedKmerIndex = -1;
        
        Iterator<Kmer> itr = leftKmers.iterator();
        for (int i=0; i<=maxLeftKmersIndex; ++i) {
            Kmer leftKmer = itr.next();
            if (rightReadPartnerIndex >=0 && rightReadPartnerIndex < numRightKmers) {
                if (graph.lookupReadKmerPair(leftKmer, rightKmers.get(rightReadPartnerIndex))) {
                    ++supportingReadKmerPairs;
                    lastSupportedKmerIndex = i;
                }
            }
            
            if (rightFragPartnerIndex >=0 && rightFragPartnerIndex < numRightKmers) {
                if (graph.lookupKmerPair(leftKmer, rightKmers.get(rightFragPartnerIndex))) {
                    ++supportingFragKmerPairs;
                    lastSupportedKmerIndex = i;
                }
            }
            
            ++rightReadPartnerIndex;
            ++rightFragPartnerIndex;
            
            if (rightReadPartnerIndex >= numRightKmers && rightFragPartnerIndex >= numRightKmers) {
                break;
            }
        }
        
        result[0] = supportingReadKmerPairs;
        result[1] = supportingFragKmerPairs;
        result[2] = lastSupportedKmerIndex;
    }

/*    
    private static ArrayDeque<Kmer> extendRightSE(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            int maxTipLen,
                                            int bound) {
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int readPairedKmersDist = graph.getReadPairedKmerDistance();
        
        final int numKmers = kmers.size();
        final float pathMinCov = getMinimumKmerCoverage(kmers, Math.max(numKmers - readPairedKmersDist, 0), numKmers);
        
        ArrayDeque<Kmer> candidates = kmers.get(numKmers-1).getSuccessors(k, numHash, graph);
        
        float bestScore = 0;
        float bestCov = 0;
        ArrayDeque<Kmer> bestExtension = null;
                
        for (Kmer candidate : candidates) {
            ArrayDeque<Kmer> e = naiveExtendRight(candidate, graph, maxTipLen, bound);
            e.addFirst(candidate);
            int readPartnerKmerIndex = numKmers - readPairedKmersDist;
            
            Iterator<Kmer> itr = e.iterator();
            
            int supportingReadKmerPairs = 0;
            int lastReadPartneredKmerIndex = -1;
                        
            int endIndex = Math.min(e.size(), readPairedKmersDist);
            for (int i=0; i<endIndex; ++i) {
                Kmer eKmer = itr.next();
                
                if (i < readPairedKmersDist && readPartnerKmerIndex >= 0) {
                    Kmer partner = kmers.get(readPartnerKmerIndex);
                    if (graph.lookupReadKmerPair(partner, eKmer)) {
                        ++supportingReadKmerPairs;
                        lastReadPartneredKmerIndex = i;
                    }
                }
                
                ++readPartnerKmerIndex;
            }
            
            if (lastReadPartneredKmerIndex >= 0) {
                float cov = getMedianKmerCoverage(e);
                float score = Math.min(pathMinCov, cov) * supportingReadKmerPairs / (lastReadPartneredKmerIndex+1);
                if (score > bestScore || (score == bestScore && cov > bestCov)) {
                    bestScore = score;
                    bestCov = cov;
                    bestExtension = e;
                    
                    itr = e.descendingIterator();
                    for (int i=e.size()-1; i>lastReadPartneredKmerIndex; --i) {
                        itr.next();
                        itr.remove();
                    }
                }
            }
        }
        
        return bestExtension;
    }
    
    private static ArrayDeque<Kmer> extendLeftSE(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLen,
                                            int bound) {
        // `kmers` is reversed
        
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int readPairedKmersDist = graph.getReadPairedKmerDistance();
        
        final int numKmers = kmers.size();
        final float pathMinCov = getMinimumKmerCoverage(kmers, Math.max(numKmers - readPairedKmersDist, 0), numKmers);
        
        ArrayDeque<Kmer> candidates = kmers.get(numKmers-1).getPredecessors(k, numHash, graph);
        
        float bestScore = 0;
        float bestCov = 0;
        ArrayDeque<Kmer> bestExtension = null;
        
        for (Kmer candidate : candidates) {
            ArrayDeque<Kmer> e = naiveExtendLeft(candidate, graph, maxTipLen, bound);
            e.addFirst(candidate);
            int readPartnerKmerIndex = numKmers - readPairedKmersDist;
            
            Iterator<Kmer> itr = e.iterator();
            
            int supportingReadKmerPairs = 0;
            int lastReadPartneredKmerIndex = -1;
            
            int endIndex = Math.min(e.size(), readPairedKmersDist);
            for (int i=0; i<endIndex; ++i) {
                Kmer eKmer = itr.next();
                
                if (i < readPairedKmersDist && readPartnerKmerIndex >= 0) {
                    Kmer partner = kmers.get(readPartnerKmerIndex);
                    if (graph.lookupReadKmerPair(eKmer, partner)) {
                        ++supportingReadKmerPairs;
                        lastReadPartneredKmerIndex = i;
                    }
                }
                
                ++readPartnerKmerIndex;
            }
            
            if (lastReadPartneredKmerIndex >= 0) {
                float cov = getMedianKmerCoverage(e);
                float score = Math.min(pathMinCov, cov) * supportingReadKmerPairs / (lastReadPartneredKmerIndex+1);
                if (score > bestScore || (score == bestScore && cov > bestCov)) {
                    bestScore = score;
                    bestCov = cov;
                    bestExtension = e;
                    
                    itr = e.descendingIterator();
                    for (int i=e.size()-1; i>lastReadPartneredKmerIndex; --i) {
                        itr.next();
                        itr.remove();
                    }
                }
            }
        }
        
        return bestExtension;
    }
*/
    
    private static ArrayDeque<Kmer> extendRightSE(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            int maxTipLen,
                                            float minKmerCov) {
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int readPairedKmersDist = graph.getReadPairedKmerDistance();
        final int numKmers = kmers.size();
        int maxExtensionLength = readPairedKmersDist - 2; // -1 for candidate k-mer; -1 for partner kmer on current sequence
        
        ArrayDeque<Kmer> candidates = kmers.get(numKmers-1).getSuccessors(k, numHash, graph);
        
        if (candidates.size() == 1) {
            Kmer c = candidates.peek();
            ArrayDeque<Kmer> e = naiveExtendRightNoBackChecks(c, graph, maxTipLen, maxExtensionLength, minKmerCov);
            e.addFirst(c);
            return e;
        }
    
        final float pathMinCov = getMinimumKmerCoverage(kmers, Math.max(numKmers - readPairedKmersDist, 0), numKmers);        
        float bestScore = 0;
        float bestCov = 0;
        ArrayDeque<Kmer> bestExtension = null;
        int[] result = new int[2];
        
        for (Kmer candidate : candidates) {
            ArrayDeque<Kmer> e = naiveExtendRightNoBackChecks(candidate, graph, maxTipLen, maxExtensionLength, minKmerCov);
            e.addFirst(candidate);
            
            countKmerPairsSE(kmers, e, 0, graph, result);
            
            int lastPartneredKmerIndex = result[1];
            
            if (lastPartneredKmerIndex >= 0 && result[0] > 0) {
                float cov = getMedianKmerCoverage(e);
                float score = Math.min(pathMinCov, cov) * result[0] / (lastPartneredKmerIndex+1);
//                System.out.println(score + ": " + graph.assemble(e));
                if (score > bestScore || (score == bestScore && cov > bestCov)) {
                    bestScore = score;
                    bestCov = cov;
                    bestExtension = e;
                    
                    Iterator<Kmer> itr = e.descendingIterator();
                    for (int i=e.size()-1; i>lastPartneredKmerIndex; --i) {
                        itr.next();
                        itr.remove();
                    }
                }
            }
            else {
                int gap = e.size();                
                
                if (gap >= readPairedKmersDist-1 && result[0] == 0) {
                    continue;
                }
                
                // not enough supporting paired k-mers in first extension
                
                ArrayDeque<Kmer> nextCandidates = e.getLast().getSuccessors(k, numHash, graph);
                
                for (Kmer nextCandidate : nextCandidates) {
                    ArrayDeque<Kmer> ne = naiveExtendRightNoBackChecks(nextCandidate, graph, maxTipLen, readPairedKmersDist-gap, minKmerCov);
                    ne.addFirst(nextCandidate);
                    Iterator<Kmer> itr = e.descendingIterator();
                    while (itr.hasNext()) {
                        ne.addFirst(itr.next());
                    }
                    
                    countKmerPairsSE(kmers, ne, 0, graph, result);
                    
                    lastPartneredKmerIndex = result[1];
                    if (lastPartneredKmerIndex >= 0 && result[0] > 0) {
                        float cov = getMedianKmerCoverage(ne);
                        float score = Math.min(pathMinCov, cov) * result[0] / (lastPartneredKmerIndex+1);

                        if (score > bestScore || (score == bestScore && cov > bestCov)) {
                            bestScore = score;
                            bestCov = cov;
                            bestExtension = ne;

                            itr = ne.descendingIterator();
                            for (int i=ne.size()-1; i>lastPartneredKmerIndex; --i) {
                                itr.next();
                                itr.remove();
                            }
                        }
                    }
                }
            }
        }
        
        return bestExtension;
    }
        
    private static ArrayDeque<Kmer> extendLeftSE(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLen,
                                            float minKmerCov) {
        // `kmers` is reversed
        
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int readPairedKmersDist = graph.getReadPairedKmerDistance();
        final int numKmers = kmers.size();
        int maxExtensionLength = readPairedKmersDist - 2; // -1 for candidate k-mer; -1 for partner kmer on current sequence
        
        ArrayDeque<Kmer> candidates = kmers.get(numKmers-1).getPredecessors(k, numHash, graph);
        
        if (candidates.size() == 1) {
            Kmer c = candidates.peek();
            ArrayDeque<Kmer> e = naiveExtendLeftNoBackChecks(c, graph, maxTipLen, maxExtensionLength, minKmerCov);
            e.addFirst(c);
            return e;
        }
        
        final float pathMinCov = getMinimumKmerCoverage(kmers, Math.max(numKmers - readPairedKmersDist, 0), numKmers);
        float bestScore = 0;
        float bestCov = 0;
        ArrayDeque<Kmer> bestExtension = null;
        int[] result = new int[2];
        
        for (Kmer candidate : candidates) {
            ArrayDeque<Kmer> e = naiveExtendLeftNoBackChecks(candidate, graph, maxTipLen, maxExtensionLength, minKmerCov);
            e.addFirst(candidate);
            
            countKmerPairsReversedSE(e, kmers, 0, graph, result);
            
            int lastPartneredKmerIndex = result[1];
            
            if (lastPartneredKmerIndex >= 0 && result[0] > 0) {
                float cov = getMedianKmerCoverage(e);
                float score = Math.min(pathMinCov, cov) * result[0] / (lastPartneredKmerIndex+1);
                if (score > bestScore || (score == bestScore && cov > bestCov)) {
                    bestScore = score;
                    bestCov = cov;
                    bestExtension = e;
                    
                    Iterator<Kmer> itr = e.descendingIterator();
                    for (int i=e.size()-1; i>lastPartneredKmerIndex; --i) {
                        itr.next();
                        itr.remove();
                    }
                }
            }
            else {
                int gap = e.size();
                
                if (gap >= readPairedKmersDist-1 && result[0] == 0) {
                    continue;
                }

                // not enough supporting paired k-mers in first extension
                
                ArrayDeque<Kmer> nextCandidates = e.getLast().getPredecessors(k, numHash, graph);
                
                for (Kmer nextCandidate : nextCandidates) {
                    ArrayDeque<Kmer> ne = naiveExtendLeftNoBackChecks(nextCandidate, graph, maxTipLen, readPairedKmersDist-gap, minKmerCov);
                    ne.addFirst(nextCandidate);
                    Iterator<Kmer> itr = e.descendingIterator();
                    while (itr.hasNext()) {
                        ne.addFirst(itr.next());
                    }
                                        
                    countKmerPairsReversedSE(ne, kmers, 0, graph, result);
                    
                    lastPartneredKmerIndex = result[1];
                    if (lastPartneredKmerIndex >= 0 && result[0] > 0) {
                        float cov = getMedianKmerCoverage(ne);
                        float score = Math.min(pathMinCov, cov) * result[0] / (lastPartneredKmerIndex+1);
                        if (score > bestScore || (score == bestScore && cov > bestCov)) {
                            bestScore = score;
                            bestCov = cov;
                            bestExtension = ne;

                            itr = ne.descendingIterator();
                            for (int i=ne.size()-1; i>lastPartneredKmerIndex; --i) {
                                itr.next();
                                itr.remove();
                            }
                        }
                    }
                }
            }
        }
        
        return bestExtension;
    }
    
    private static ArrayDeque<Kmer> extendRightPE(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            int maxTipLen,
                                            float minKmerCov) {
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int readPairedKmersDist = graph.getReadPairedKmerDistance();
        final int fragPairedKmersDist = graph.getFragPairedKmerDistance();        
        final int numKmers = kmers.size();
        int maxExtensionLength = fragPairedKmersDist - 2; // -1 for candidate k-mer; -1 for partner kmer on current sequence
        
        ArrayDeque<Kmer> candidates = kmers.get(numKmers-1).getSuccessors(k, numHash, graph);
        
        if (candidates.size() == 1) {
            Kmer c = candidates.peek();
            ArrayDeque<Kmer> e = naiveExtendRightNoBackChecks(c, graph, maxTipLen, maxExtensionLength, minKmerCov);
            e.addFirst(c);
            return e;
        }
    
        for (int i=numKmers-1; i>=0; --i) {
            if (graph.isRepeatKmer(kmers.get(i))) {
                --maxExtensionLength;
            }
            else {
               break;
            }
        }
        
        final float pathMinCov = getMinimumKmerCoverage(kmers, Math.max(numKmers - fragPairedKmersDist, 0), numKmers);        
        float bestScore = 0;
        float bestCov = 0;
        ArrayDeque<Kmer> bestExtension = null;
        int[] result = new int[3];
        
        for (Kmer candidate : candidates) {
            ArrayDeque<Kmer> e = naiveExtendRightNoBackChecks(candidate, graph, maxTipLen, maxExtensionLength, minKmerCov);
            e.addFirst(candidate);
            
            countKmerPairsPE(kmers, e, 0, graph, result);
            
            int lastPartneredKmerIndex = result[2];
            
            if (lastPartneredKmerIndex >= 0 && result[0] > 0 && result[1] > 0) {
                float cov = getMedianKmerCoverage(e);
                float score = Math.min(pathMinCov, cov) * (result[0] + result[1]) / (lastPartneredKmerIndex+1);
//                System.out.println(score + ": " + graph.assemble(e));
                if (score > bestScore || (score == bestScore && cov > bestCov)) {
                    bestScore = score;
                    bestCov = cov;
                    bestExtension = e;
                    
                    Iterator<Kmer> itr = e.descendingIterator();
                    for (int i=e.size()-1; i>lastPartneredKmerIndex; --i) {
                        itr.next();
                        itr.remove();
                    }
                }
            }
            else {
                int gap = e.size();                
                
                if ((gap >= readPairedKmersDist-1 && result[0] == 0) ||
                        (gap >= fragPairedKmersDist-1 && result[1] == 0)) {
                    continue;
                }
                
                // not enough supporting paired k-mers in first extension
                
                ArrayDeque<Kmer> nextCandidates = e.getLast().getSuccessors(k, numHash, graph);
                
                for (Kmer nextCandidate : nextCandidates) {
                    ArrayDeque<Kmer> ne = naiveExtendRightNoBackChecks(nextCandidate, graph, maxTipLen, maxExtensionLength-gap, minKmerCov);
                    ne.addFirst(nextCandidate);
                    Iterator<Kmer> itr = e.descendingIterator();
                    while (itr.hasNext()) {
                        ne.addFirst(itr.next());
                    }
                    
                    countKmerPairsPE(kmers, ne, 0, graph, result);
                    
                    lastPartneredKmerIndex = result[2];
                    if (lastPartneredKmerIndex >= 0 && result[0] > 0 && result[1] > 0) {
                        float cov = getMedianKmerCoverage(ne);
                        float score = Math.min(pathMinCov, cov) * (result[0] + result[1]) / (lastPartneredKmerIndex+1);

                        if (score > bestScore || (score == bestScore && cov > bestCov)) {
                            bestScore = score;
                            bestCov = cov;
                            bestExtension = ne;

                            itr = ne.descendingIterator();
                            for (int i=ne.size()-1; i>lastPartneredKmerIndex; --i) {
                                itr.next();
                                itr.remove();
                            }
                        }
                    }
                }
            }
        }
        
        return bestExtension;
    }
        
    private static ArrayDeque<Kmer> extendLeftPE(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLen,
                                            float minKmerCov) {
        // `kmers` is reversed
        
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int readPairedKmersDist = graph.getReadPairedKmerDistance();
        final int fragPairedKmersDist = graph.getFragPairedKmerDistance();
        final int numKmers = kmers.size();
        int maxExtensionLength = fragPairedKmersDist - 2; // -1 for candidate k-mer; -1 for partner kmer on current sequence
        
        ArrayDeque<Kmer> candidates = kmers.get(numKmers-1).getPredecessors(k, numHash, graph);
        
        if (candidates.size() == 1) {
            Kmer c = candidates.peek();
            ArrayDeque<Kmer> e = naiveExtendLeftNoBackChecks(c, graph, maxTipLen, maxExtensionLength, minKmerCov);
            e.addFirst(c);
            return e;
        }
        
        for (int i=numKmers-1; i>=0; --i) {
            if (graph.isRepeatKmer(kmers.get(i))) {
                --maxExtensionLength;
            }
            else {
               break;
            }
        }
        
        final float pathMinCov = getMinimumKmerCoverage(kmers, Math.max(numKmers - fragPairedKmersDist, 0), numKmers);
        float bestScore = 0;
        float bestCov = 0;
        ArrayDeque<Kmer> bestExtension = null;
        int[] result = new int[3];
        
        for (Kmer candidate : candidates) {
            ArrayDeque<Kmer> e = naiveExtendLeftNoBackChecks(candidate, graph, maxTipLen, maxExtensionLength, minKmerCov);
            e.addFirst(candidate);
            
            countKmerPairsReversedPE(e, kmers, 0, graph, result);
            
            int lastPartneredKmerIndex = result[2];
            
            if (lastPartneredKmerIndex >= 0 && result[0] > 0 && result[1] > 0) {
                float cov = getMedianKmerCoverage(e);
                float score = Math.min(pathMinCov, cov) * (result[0] + result[1]) / (lastPartneredKmerIndex+1);
                if (score > bestScore || (score == bestScore && cov > bestCov)) {
                    bestScore = score;
                    bestCov = cov;
                    bestExtension = e;
                    
                    Iterator<Kmer> itr = e.descendingIterator();
                    for (int i=e.size()-1; i>lastPartneredKmerIndex; --i) {
                        itr.next();
                        itr.remove();
                    }
                }
            }
            else {
                int gap = e.size();
                
                if ((gap >= readPairedKmersDist-1 && result[0] == 0) ||
                        (gap >= fragPairedKmersDist-1 && result[1] == 0)) {
                    continue;
                }

                // not enough supporting paired k-mers in first extension
                
                ArrayDeque<Kmer> nextCandidates = e.getLast().getPredecessors(k, numHash, graph);
                
                for (Kmer nextCandidate : nextCandidates) {
                    ArrayDeque<Kmer> ne = naiveExtendLeftNoBackChecks(nextCandidate, graph, maxTipLen, maxExtensionLength-gap, minKmerCov);
                    ne.addFirst(nextCandidate);
                    Iterator<Kmer> itr = e.descendingIterator();
                    while (itr.hasNext()) {
                        ne.addFirst(itr.next());
                    }
                                        
                    countKmerPairsReversedPE(ne, kmers, 0, graph, result);
                    
                    lastPartneredKmerIndex = result[2];
                    if (lastPartneredKmerIndex >= 0 && result[0] > 0 && result[1] > 0) {
                        float cov = getMedianKmerCoverage(ne);
                        float score = Math.min(pathMinCov, cov) * (result[0] + result[1]) / (lastPartneredKmerIndex+1);
                        if (score > bestScore || (score == bestScore && cov > bestCov)) {
                            bestScore = score;
                            bestCov = cov;
                            bestExtension = ne;

                            itr = ne.descendingIterator();
                            for (int i=ne.size()-1; i>lastPartneredKmerIndex; --i) {
                                itr.next();
                                itr.remove();
                            }
                        }
                    }
                }
            }
        }
        
        return bestExtension;
    }

    private static boolean hasDuplicatedKmerPair(ArrayList<Kmer> kmers, Kmer cursor, int d, int mateIndex) {
        if (mateIndex >= 0) {
            Kmer mate = kmers.get(mateIndex);

            int cursor2Index = kmers.lastIndexOf(cursor);
            if (cursor2Index >= 0) {
                int mate2Index = cursor2Index - d;

                if (mate2Index >= 0) {
                    if (mate.equals(kmers.get(mate2Index))) {
                        return true;
                    }

                    int cursor1Index = kmers.indexOf(cursor);
                    if (cursor1Index != cursor2Index) {
                        int mate1Index = cursor1Index - d;

                        if (mate1Index >= 0 && mate.equals(kmers.get(mate1Index))) {
                            return true;
                        }

                        for (int i=cursor1Index+1; i<cursor2Index; ++i) {
                            Kmer candidate = kmers.get(i);
                            if (cursor.equals(candidate)) {
                                int i2 = i - d;
                                if (i2 >= 0 && mate.equals(kmers.get(i2))) {
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        return false;
    }
    
    public static int[] extendSE(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            int maxTipLength,
                                            float minKmerCov) {
        final int d = graph.getReadPairedKmerDistance();
        final float multiplier = 0.1f;
        
        int origLen = kmers.size();
        
        HashSet<Kmer> usedKmers = new HashSet<>(kmers);

        Collections.reverse(kmers);
        
        while (true) {
            float covThreshold = getMinimumKmerCoverage(kmers, Math.max(0, kmers.size()-d), kmers.size());
            
            ArrayDeque<Kmer> e =  null;
                
            while (true) {
                covThreshold = Math.max(minKmerCov, covThreshold * multiplier);
                
                e =  extendLeftSE(kmers, graph, maxTipLength, covThreshold);
                
                if ((e != null && !e.isEmpty()) || covThreshold == minKmerCov) {
                    break;
                }
            }

            if (e == null || e.isEmpty()) {
                break;
            }
            
            Iterator<Kmer> itr = e.iterator();
            boolean used = true;
            
            while (itr.hasNext() && used) {
                Kmer kmer = itr.next();
                
                if (!usedKmers.contains(kmer)) {
                    used = false;
                }
            }
            
            int endIndex = Math.max(0, kmers.size()-d+e.size());
            for (int i=kmers.size()-1; i>=endIndex && used; --i) {
                Kmer kmer = kmers.get(i);
                
                if (!usedKmers.contains(kmer)) {
                    used = false;
                }
            }
            
            if (used && (kmers.size() < d || hasDuplicatedKmerPair(kmers, e.getLast(), d, kmers.size()-1-d+e.size()))) {
                break;
            }
            
            kmers.addAll(e);
            usedKmers.addAll(e);
        }
        
        int leftExtLen = kmers.size() - origLen;
        
        Collections.reverse(kmers);
        
        while (true) {
            float covThreshold = getMinimumKmerCoverage(kmers, Math.max(0, kmers.size()-d), kmers.size());
            
            ArrayDeque<Kmer> e =  null;
                
            while (true) {
                covThreshold = Math.max(minKmerCov, covThreshold * multiplier);
                
                e =  extendRightSE(kmers, graph, maxTipLength, covThreshold);
                
                if ((e != null && !e.isEmpty()) || covThreshold == minKmerCov) {
                    break;
                }
            }

            if (e == null || e.isEmpty()) {
                break;
            }
            
            Iterator<Kmer> itr = e.iterator();
            boolean used = true;
            while (itr.hasNext() && used) {
                Kmer kmer = itr.next();
                
                if (!usedKmers.contains(kmer)) {
                    used = false;
                }
            }
            
            int endIndex = Math.max(0, kmers.size()-d+e.size());
            for (int i=kmers.size()-1; i>=endIndex && used; --i) {
                Kmer kmer = kmers.get(i);
                
                if (!usedKmers.contains(kmer)) {
                    used = false;
                }
            }
            
            if (used && (kmers.size() < d || hasDuplicatedKmerPair(kmers, e.getLast(), d, kmers.size()-1-d+e.size()))) {
                break;
            }
            
            kmers.addAll(e);
            usedKmers.addAll(e);
        }
                
        return new int[]{leftExtLen, leftExtLen + origLen};
    }
    
    public static int[] extendPE(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph,
                                            int maxTipLength,
                                            float minKmerCov) {
        final int d = graph.getFragPairedKmerDistance();
        final float multiplier = 0.1f;
        
        int origLen = kmers.size();
        
        HashSet<Kmer> usedKmers = new HashSet<>(kmers);

        Collections.reverse(kmers);
        
        while (true) {
            float covThreshold = getMinimumKmerCoverage(kmers, Math.max(0, kmers.size()-d), kmers.size());
            
            ArrayDeque<Kmer> e =  null;
                
            while (true) {
                covThreshold = Math.max(minKmerCov, covThreshold * multiplier);
                
                e =  extendLeftPE(kmers, graph, maxTipLength, covThreshold);
                
                if ((e != null && !e.isEmpty()) || covThreshold == minKmerCov) {
                    break;
                }
            }

            if (e == null || e.isEmpty()) {
                break;
            }
            
            Iterator<Kmer> itr = e.iterator();
            boolean used = true;
            
            while (itr.hasNext() && used) {
                Kmer kmer = itr.next();
                
                if (!usedKmers.contains(kmer)) {
                    used = false;
                }
            }
            
            int endIndex = Math.max(0, kmers.size()-d+e.size());
            for (int i=kmers.size()-1; i>=endIndex && used; --i) {
                Kmer kmer = kmers.get(i);
                
                if (!usedKmers.contains(kmer)) {
                    used = false;
                }
            }
            
            if (used && (kmers.size() < d || hasDuplicatedKmerPair(kmers, e.getLast(), d, kmers.size()-1-d+e.size()))) {
                break;
            }
            
            kmers.addAll(e);
            usedKmers.addAll(e);
        }
        
        int leftExtLen = kmers.size() - origLen;
        
        Collections.reverse(kmers);
        
        while (true) {
            float covThreshold = getMinimumKmerCoverage(kmers, Math.max(0, kmers.size()-d), kmers.size());
            
            ArrayDeque<Kmer> e =  null;
                
            while (true) {
                covThreshold = Math.max(minKmerCov, covThreshold * multiplier);
                
                e =  extendRightPE(kmers, graph, maxTipLength, covThreshold);
                
                if ((e != null && !e.isEmpty()) || covThreshold == minKmerCov) {
                    break;
                }
            }

            if (e == null || e.isEmpty()) {
                break;
            }
            
            Iterator<Kmer> itr = e.iterator();
            boolean used = true;
            while (itr.hasNext() && used) {
                Kmer kmer = itr.next();
                
                if (!usedKmers.contains(kmer)) {
                    used = false;
                }
            }
            
            int endIndex = Math.max(0, kmers.size()-d+e.size());
            for (int i=kmers.size()-1; i>=endIndex && used; --i) {
                Kmer kmer = kmers.get(i);
                
                if (!usedKmers.contains(kmer)) {
                    used = false;
                }
            }
            
            if (used && (kmers.size() < d || hasDuplicatedKmerPair(kmers, e.getLast(), d, kmers.size()-1-d+e.size()))) {
                break;
            }
            
            kmers.addAll(e);
            usedKmers.addAll(e);
        }
                
        return new int[]{leftExtLen, leftExtLen + origLen};
    }
    
    public static boolean extendRightWithPairedKmersDFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            BloomFilter assembledKmersBloomFilter,
                                            HashSet<Kmer> usedKmers,
                                            float maxCovGradient) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        final int distance = graph.getFragPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int maxDepth = maxRightPartnerSearchDepth2(kmers, graph, distance, assembledKmersBloomFilter, minNumPairs);
        
        // data structure to store visited kmers at defined depth
        HashMap<Kmer, ArrayDeque<Integer>> visitedKmers = new HashMap<>();
                
        int numKmers = kmers.size();
        int depth = 0;
        int partnerIndex = numKmers - distance + depth;
        
//        float minEdgeCoverage = getMinimumKmerCoverage(kmers, Math.max(0, numKmers-distance), numKmers-1);
//        if (minEdgeCoverage >= k && maxDepth > k) {
//            maxDepth = Math.min(distance/2, maxDepth);
//        }
        
//        maxDepth = Math.min(maxDepth, k);
        
//        float minCovThreshold = (float) Math.floor(minEdgeCoverage * maxCovGradient);
        
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
//        branchesStack.add(getSuccessorsRanked(kmers.get(numKmers-1), graph, lookahead, minCovThreshold));
        branchesStack.add(getSuccessorsRanked(kmers.get(numKmers-1), graph, lookahead));
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        HashSet<Kmer> extensionKmers = new HashSet<>();
        
        int maxPartnerIndex = numKmers - 1 - minNumPairs;
        while (!branchesStack.isEmpty()) {
            LinkedList<Kmer> branches = branchesStack.getLast();
            
            if (branches.isEmpty()) {
                Kmer cursor = extension.pollLast();
                
                if (cursor != null) {
                    extensionKmers.remove(cursor);
                }
                
                branchesStack.removeLast();
                --depth;
                --partnerIndex;
            }
            else {
                Kmer cursor = branches.pop();
                
//                if (!graph.isLowComplexity(cursor)) {
                    if (partnerIndex >=0 &&
                            partnerIndex <= maxPartnerIndex &&
    //                        partnerIndex+minNumPairs < kmers.size() && 
                            hasPairedRightKmers(cursor, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {

                        if (usedKmers.contains(cursor)) {
                            // check if all extension kmers has been used in this sequence already
                            boolean allExtensionKmersUsed = true;
                            for (Kmer kmer : extension) {
                                if (!usedKmers.contains(kmer)) {
                                    allExtensionKmersUsed = false;
                                    break;
                                }
                            }
                            
                            if (allExtensionKmersUsed) {
                                // check whether this kmer pair has been used in this sequence already
                                Kmer partner = kmers.get(partnerIndex);
                                for (int i = kmers.size()-1; i >= distance; --i) {
                                    if (cursor.equals(kmers.get(i)) &&
                                            partner.equals(kmers.get(i-distance))) {
                                        return false;
                                    }
                                }
                            }
                        }

                        if (assembledKmersBloomFilter.lookup(cursor.getHash())) {
                            boolean assembled = greedyExtendRight(graph, cursor, lookahead, minNumPairs, assembledKmersBloomFilter) != null;

                            if (assembled) {
                                for (Kmer kmer : extension) {
                                    if (!assembledKmersBloomFilter.lookup(kmer.getHash())) {
                                            assembled = false;
                                            break;
                                    }
                                }

                                if (assembled) {
                                    for (int i=partnerIndex; i<numKmers; ++i) {
                                        if (!assembledKmersBloomFilter.lookup(kmers.get(i).getHash())) {
                                            assembled = false;
                                            break;
                                        }
                                    }
                                }
                            }

                            if (assembled) {
                                // region already assembled, do not extend further
                                return false;
                            }
                        }

                        kmers.addAll(extension);
                        kmers.add(cursor);

                        usedKmers.addAll(extension);
                        usedKmers.add(cursor);

                        return true;
                    }
                    else if (depth < maxDepth &&
                            (depth == 0 || !extensionKmers.contains(cursor))) {

//                        float minExtensionCoverage = extension.isEmpty() ? minEdgeCoverage : Math.min(minEdgeCoverage, getMinimumKmerCoverage(extension));
//                        minCovThreshold = (float) Math.floor(minExtensionCoverage * maxCovGradient);

                        if (cursor.hasAtLeastXPredecessors(k, numHash, graph, 2)) {
                            // only consider kmers that may be visited from an alternative branch upstream

                            ArrayDeque<Integer> visitedDepths = visitedKmers.get(cursor);
                            if (visitedDepths == null) {
                                visitedDepths = new ArrayDeque<>();
                                visitedDepths.add(depth);

                                visitedKmers.put(cursor, visitedDepths);

//                                branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead, minCovThreshold));
                                branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursor);
                                ++depth;
                                ++partnerIndex;
                            }
                            else {
                                boolean visited = false;
                                for (Integer d : visitedDepths) {
                                    if (Math.abs(d - depth) <= maxIndelSize) {
                                        visited = true;
                                        break;
                                    }
                                }

                                if (!visited) {
                                    visitedDepths.add(depth);

//                                    branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead, minCovThreshold));
                                    branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                                    extension.add(cursor);
                                    extensionKmers.add(cursor);
                                    ++depth;
                                    ++partnerIndex;
                                }
                            }
                        }
                        else {
//                            branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead, minCovThreshold));
                            branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursor);
                            ++depth;
                            ++partnerIndex;
                        }
                    }
//                }
            }
        }
        
        return false;
    }
    
    public static boolean extendLeftWithPairedKmersDFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            BloomFilter assembledKmersBloomFilter,
                                            HashSet<Kmer> usedKmers,
                                            float maxCovGradient) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        final int distance = graph.getFragPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int maxDepth = maxLeftPartnerSearchDepth2(kmers, graph, distance, assembledKmersBloomFilter, minNumPairs);
        
        // data structure to store visited kmers at defined depth
        HashMap<Kmer, ArrayDeque<Integer>> visitedKmers = new HashMap<>();
        
        int numKmers = kmers.size();
        int depth = 0;
        int partnerIndex = numKmers - distance + depth;
        
//        float minEdgeCoverage = getMinimumKmerCoverage(kmers, Math.max(0, numKmers-distance), numKmers-1);
//        if (minEdgeCoverage >= k && maxDepth > k) {
//            maxDepth = Math.min(distance/2, maxDepth);
//        }
//        maxDepth = Math.min(maxDepth, k);

//        float minCovThreshold = (float) Math.floor(minEdgeCoverage * maxCovGradient);
        
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
//        branchesStack.add(getPredecessorsRanked(kmers.get(numKmers-1), graph, lookahead, minCovThreshold));
        branchesStack.add(getPredecessorsRanked(kmers.get(numKmers-1), graph, lookahead));
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        HashSet<Kmer> extensionKmers = new HashSet<>();
        
        int maxPartnerIndex = numKmers - 1 - minNumPairs;
        while (!branchesStack.isEmpty()) {
            LinkedList<Kmer> branches = branchesStack.getLast();
            
            if (branches.isEmpty()) {
                Kmer cursor = extension.pollLast();
                if (cursor != null) {
                    extensionKmers.remove(cursor);
                }
                
                branchesStack.removeLast();
                --depth;
                --partnerIndex;
            }
            else {
                Kmer cursor = branches.pop();
                
//                if (!graph.isLowComplexity(cursor)) {
                    if (partnerIndex >=0 &&
                            partnerIndex <= maxPartnerIndex &&
    //                        partnerIndex+minNumPairs < kmers.size() && 
                            hasPairedLeftKmers(cursor, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {

                        if (usedKmers.contains(cursor)) {
                            // check if all extension kmers has been used in this sequence already
                            boolean allExtensionKmersUsed = true;
                            for (Kmer kmer : extension) {
                                if (!usedKmers.contains(kmer)) {
                                    allExtensionKmersUsed = false;
                                    break;
                                }
                            }
                            
                            if (allExtensionKmersUsed) {
                                // check whether this kmer pair has been used in this sequence already (ie. a loop)
                                Kmer partner = kmers.get(partnerIndex);
                                for (int i = kmers.size()-1; i >= distance; --i) {
                                    if (cursor.equals(kmers.get(i)) &&
                                            partner.equals(kmers.get(i-distance))) {
                                        return false;
                                    }
                                }
                            }
                        }

                        if (assembledKmersBloomFilter.lookup(cursor.getHash())) {
                            boolean assembled = greedyExtendLeft(graph, cursor, lookahead, minNumPairs, assembledKmersBloomFilter) != null;

                            if (assembled) {
                                // check if all extension kmers has been used
                                for (Kmer kmer : extension) {
                                    if (!assembledKmersBloomFilter.lookup(kmer.getHash())) {
                                        assembled = false;
                                        break;
                                    }
                                }

                                if (assembled) {
                                    for (int i=partnerIndex; i<numKmers; ++i) {
                                        if (!assembledKmersBloomFilter.lookup(kmers.get(i).getHash())) {
    //                                        if (distanceInversePI < ++numNotAssembled) {
                                                assembled = false;
                                                break;
    //                                        }
                                        }
                                    }
                                }
                            }

                            if (assembled) {
                                // region already assembled, do not extend further
                                return false;
                            }
                        }

                        kmers.addAll(extension);
                        kmers.add(cursor);

                        usedKmers.addAll(extension);
                        usedKmers.add(cursor);

                        return true;
                    }
                    else if (depth < maxDepth &&
                            (depth == 0 || !extensionKmers.contains(cursor))) {

//                        float minExtensionCoverage = extension.isEmpty() ? minEdgeCoverage : Math.min(minEdgeCoverage, getMinimumKmerCoverage(extension));
//                        minCovThreshold = (float) Math.floor(minExtensionCoverage * maxCovGradient);

                        if (cursor.hasAtLeastXSuccessors(k, numHash, graph, 2)) {
                            // only consider kmers that may be visited from an alternative branch upstream

                            ArrayDeque<Integer> visitedDepths = visitedKmers.get(cursor);
                            if (visitedDepths == null) {
                                visitedDepths = new ArrayDeque<>();
                                visitedDepths.add(depth);

                                visitedKmers.put(cursor, visitedDepths);

//                                branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead, minCovThreshold));
                                branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursor);
                                ++depth;
                                ++partnerIndex;
                            }
                            else {
                                boolean visited = false;
                                for (Integer d : visitedDepths) {
                                    if (Math.abs(d - depth) <= maxIndelSize) {
                                        visited = true;
                                        break;
                                    }
                                }

                                if (!visited) {
                                    visitedDepths.add(depth);

//                                    branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead, minCovThreshold));
                                    branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                                    extension.add(cursor);
                                    extensionKmers.add(cursor);
                                    ++depth;
                                    ++partnerIndex;
                                }
                            }
                        }
                        else {
//                            branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead, minCovThreshold));
                            branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursor);
                            ++depth;
                            ++partnerIndex;
                        }
                    }
//                }
            }
        }
        
        return false;
    }
    
    public static boolean extendRightWithPairedKmersDFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            HashSet<Kmer> usedKmers) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        final int distance = graph.getFragPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int maxDepth = maxRightPartnerSearchDepth2(kmers, graph, distance, minNumPairs);
        
        // data structure to store visited kmers at defined depth
        HashMap<Kmer, ArrayDeque<Integer>> visitedKmers = new HashMap<>();
                
        int numKmers = kmers.size();
        int depth = 0;
        int partnerIndex = numKmers - distance + depth;
        
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
        branchesStack.add(getSuccessorsRanked(kmers.get(numKmers-1), graph, lookahead));
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        HashSet<Kmer> extensionKmers = new HashSet<>();
        
        int maxPartnerIndex = numKmers - 1 - minNumPairs;
        while (!branchesStack.isEmpty()) {
            LinkedList<Kmer> branches = branchesStack.getLast();
            
            if (branches.isEmpty()) {
                Kmer cursor = extension.pollLast();
                
                if (cursor != null) {
                    extensionKmers.remove(cursor);
                }
                
                branchesStack.removeLast();
                --depth;
                --partnerIndex;
            }
            else {
                Kmer cursor = branches.pop();
                
                if (!graph.isLowComplexity(cursor)) {
                    if (partnerIndex >=0 &&
                            partnerIndex <= maxPartnerIndex &&
    //                        partnerIndex+minNumPairs < kmers.size() && 
                            hasPairedRightKmers(cursor, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {

                        if (usedKmers.contains(cursor)) {
                            // check whether this kmer pair has been used in this sequence already
                            Kmer partner = kmers.get(partnerIndex);
                            for (int i = kmers.size()-1; i >= distance; --i) {
                                if (cursor.equals(kmers.get(i)) &&
                                        partner.equals(kmers.get(i-distance))) {
                                    return false;
                                }
                            }
                        }

                        kmers.addAll(extension);
                        kmers.add(cursor);

                        usedKmers.addAll(extension);
                        usedKmers.add(cursor);

                        return true;
                    }
                    else if (depth < maxDepth &&
                            (depth == 0 || !extensionKmers.contains(cursor))) {

                        if (cursor.hasAtLeastXPredecessors(k, numHash, graph, 2)) {
                            // only consider kmers that may be visited from an alternative branch upstream

                            ArrayDeque<Integer> visitedDepths = visitedKmers.get(cursor);
                            if (visitedDepths == null) {
                                visitedDepths = new ArrayDeque<>();
                                visitedDepths.add(depth);

                                visitedKmers.put(cursor, visitedDepths);

                                branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursor);
                                ++depth;
                                ++partnerIndex;
                            }
                            else {
                                boolean visited = false;
                                for (Integer d : visitedDepths) {
                                    if (Math.abs(d - depth) <= maxIndelSize) {
                                        visited = true;
                                        break;
                                    }
                                }

                                if (!visited) {
                                    visitedDepths.add(depth);

                                    branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                                    extension.add(cursor);
                                    extensionKmers.add(cursor);
                                    ++depth;
                                    ++partnerIndex;
                                }
                            }
                        }
                        else {
                            branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursor);
                            ++depth;
                            ++partnerIndex;
                        }
                    }
                }
            }
        }
        
        return false;
    }
    
    public static boolean extendLeftWithPairedKmersDFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            HashSet<Kmer> usedKmers) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        final int distance = graph.getFragPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int maxDepth = maxLeftPartnerSearchDepth2(kmers, graph, distance, minNumPairs);
        
        // data structure to store visited kmers at defined depth
        HashMap<Kmer, ArrayDeque<Integer>> visitedKmers = new HashMap<>();
        
        int numKmers = kmers.size();
        int depth = 0;
        int partnerIndex = numKmers - distance + depth;
                
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
        branchesStack.add(getPredecessorsRanked(kmers.get(numKmers-1), graph, lookahead));
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        HashSet<Kmer> extensionKmers = new HashSet<>();
        
        int maxPartnerIndex = numKmers - 1 - minNumPairs;
        while (!branchesStack.isEmpty()) {
            LinkedList<Kmer> branches = branchesStack.getLast();
            
            if (branches.isEmpty()) {
                Kmer cursor = extension.pollLast();
                if (cursor != null) {
                    extensionKmers.remove(cursor);
                }
                
                branchesStack.removeLast();
                --depth;
                --partnerIndex;
            }
            else {
                Kmer cursor = branches.pop();
                
                if (!graph.isLowComplexity(cursor)) {
                    if (partnerIndex >=0 &&
                            partnerIndex <= maxPartnerIndex &&
    //                        partnerIndex+minNumPairs < kmers.size() && 
                            hasPairedLeftKmers(cursor, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {

                        if (usedKmers.contains(cursor)) {
                            // check whether this kmer pair has been used in this sequence already (ie. a loop)
                            Kmer partner = kmers.get(partnerIndex);
                            for (int i = kmers.size()-1; i >= distance; --i) {
                                if (cursor.equals(kmers.get(i)) &&
                                        partner.equals(kmers.get(i-distance))) {
                                    return false;
                                }
                            }
                        }

                        kmers.addAll(extension);
                        kmers.add(cursor);

                        usedKmers.addAll(extension);
                        usedKmers.add(cursor);

                        return true;
                    }
                    else if (depth < maxDepth &&
                            (depth == 0 || !extensionKmers.contains(cursor))) {

                        if (cursor.hasAtLeastXSuccessors(k, numHash, graph, 2)) {
                            // only consider kmers that may be visited from an alternative branch upstream

                            ArrayDeque<Integer> visitedDepths = visitedKmers.get(cursor);
                            if (visitedDepths == null) {
                                visitedDepths = new ArrayDeque<>();
                                visitedDepths.add(depth);

                                visitedKmers.put(cursor, visitedDepths);

                                branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursor);
                                ++depth;
                                ++partnerIndex;
                            }
                            else {
                                boolean visited = false;
                                for (Integer d : visitedDepths) {
                                    if (Math.abs(d - depth) <= maxIndelSize) {
                                        visited = true;
                                        break;
                                    }
                                }

                                if (!visited) {
                                    visitedDepths.add(depth);

                                    branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                                    extension.add(cursor);
                                    extensionKmers.add(cursor);
                                    ++depth;
                                    ++partnerIndex;
                                }
                            }
                        }
                        else {
                            branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursor);
                            ++depth;
                            ++partnerIndex;
                        }
                    }
                }
            }
        }
        
        return false;
    }
    
    public static boolean extendRightWithPairedKmersDFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            BloomFilter assembledKmersBloomFilter,
                                            HashSet<Kmer> usedKmers) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        final int distance = graph.getFragPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int maxDepth = maxRightPartnerSearchDepth2(kmers, graph, distance, assembledKmersBloomFilter, minNumPairs);
        
        // data structure to store visited kmers at defined depth
        HashMap<Kmer, ArrayDeque<Integer>> visitedKmers = new HashMap<>();
                
        int numKmers = kmers.size();
        int depth = 0;
        int partnerIndex = numKmers - distance + depth;
        
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
        
        branchesStack.add(getSuccessorsRanked(kmers.get(numKmers-1), graph, lookahead));
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        HashSet<Kmer> extensionKmers = new HashSet<>();
        
        Kmer cursor;
        int maxPartnerIndex = numKmers - 1 - minNumPairs;
        while (!branchesStack.isEmpty()) {
            LinkedList<Kmer> branches = branchesStack.getLast();
            
            if (branches.isEmpty()) {
                cursor = extension.pollLast();
                
                if (cursor != null) {
                    extensionKmers.remove(cursor);
                }
                
                branchesStack.removeLast();
                --depth;
                --partnerIndex;
            }
            else {
                cursor = branches.pop();
                
                if (!graph.isLowComplexity(cursor)) {
                    if (partnerIndex >=0 &&
                            partnerIndex <= maxPartnerIndex &&
    //                        partnerIndex+minNumPairs < kmers.size() && 
                            hasPairedRightKmers(cursor, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {

                        if (usedKmers.contains(cursor)) {
                            // check whether this kmer pair has been used in this sequence already
                            Kmer partner = kmers.get(partnerIndex);
                            for (int i = kmers.size()-1; i >= distance; --i) {
                                if (cursor.equals(kmers.get(i)) &&
                                        partner.equals(kmers.get(i-distance))) {
                                    return false;
                                }
                            }
                        }

                        if (assembledKmersBloomFilter.lookup(cursor.getHash())) {
                            boolean assembled = greedyExtendRight(graph, cursor, lookahead, lookahead, assembledKmersBloomFilter) != null;

                            if (assembled) {
    //                            int numNotAssembled = 0;

                                for (Kmer kmer : extension) {
                                    if (!assembledKmersBloomFilter.lookup(kmer.getHash())) {
    //                                    if (distanceInversePI < ++numNotAssembled) {
                                            assembled = false;
                                            break;
    //                                    }
                                    }
                                }

                                if (assembled) {
                                    for (int i=partnerIndex; i<numKmers; ++i) {
                                        if (!assembledKmersBloomFilter.lookup(kmers.get(i).getHash())) {
    //                                        if (distanceInversePI < ++numNotAssembled) {
                                                assembled = false;
                                                break;
    //                                        }
                                        }
                                    }
                                }
                            }

                            if (assembled) {
                                // region already assembled, do not extend further
                                return false;
                            }
                        }

                        kmers.addAll(extension);
                        kmers.add(cursor);

                        usedKmers.addAll(extension);
                        usedKmers.add(cursor);

                        return true;
                    }
                    else if (depth < maxDepth &&
                            (depth == 0 || !extensionKmers.contains(cursor))) {

                        if (cursor.hasAtLeastXPredecessors(k, numHash, graph, 2)) {
                            // only consider kmers that may be visited from an alternative branch upstream

                            ArrayDeque<Integer> visitedDepths = visitedKmers.get(cursor);
                            if (visitedDepths == null) {
                                visitedDepths = new ArrayDeque<>();
                                visitedDepths.add(depth);

                                visitedKmers.put(cursor, visitedDepths);

                                branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursor);
                                ++depth;
                                ++partnerIndex;
                            }
                            else {
                                boolean visited = false;
                                for (Integer d : visitedDepths) {
                                    if (Math.abs(d - depth) <= maxIndelSize) {
                                        visited = true;
                                        break;
                                    }
                                }

                                if (!visited) {
                                    visitedDepths.add(depth);

                                    branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                                    extension.add(cursor);
                                    extensionKmers.add(cursor);
                                    ++depth;
                                    ++partnerIndex;
                                }
                            }
                        }
                        else {
                            branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursor);
                            ++depth;
                            ++partnerIndex;
                        }
                    }
                }
            }
        }
        
        return false;
    }
    
    public static boolean extendLeftWithPairedKmersDFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            BloomFilter assembledKmersBloomFilter,
                                            HashSet<Kmer> usedKmers) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        final int distance = graph.getFragPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int maxDepth = maxLeftPartnerSearchDepth2(kmers, graph, distance, assembledKmersBloomFilter, minNumPairs);
        
        // data structure to store visited kmers at defined depth
        HashMap<Kmer, ArrayDeque<Integer>> visitedKmers = new HashMap<>();
        
        int numKmers = kmers.size();
        int depth = 0;
        int partnerIndex = numKmers - distance + depth;
        
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
        
        branchesStack.add(getPredecessorsRanked(kmers.get(numKmers-1), graph, lookahead));
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        HashSet<Kmer> extensionKmers = new HashSet<>();
        
        Kmer cursor;
        int maxPartnerIndex = numKmers - 1 - minNumPairs;
        while (!branchesStack.isEmpty()) {
            LinkedList<Kmer> branches = branchesStack.getLast();
            
            if (branches.isEmpty()) {
                cursor = extension.pollLast();
                if (cursor != null) {
                    extensionKmers.remove(cursor);
                }
                
                branchesStack.removeLast();
                --depth;
                --partnerIndex;
            }
            else {
                cursor = branches.pop();
                
                if (!graph.isLowComplexity(cursor)) {
                    if (partnerIndex >=0 &&
                            partnerIndex <= maxPartnerIndex &&
    //                        partnerIndex+minNumPairs < kmers.size() && 
                            hasPairedLeftKmers(cursor, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {

                        if (usedKmers.contains(cursor)) {
                            // check whether this kmer pair has been used in this sequence already (ie. a loop)
                            Kmer partner = kmers.get(partnerIndex);
                            for (int i = kmers.size()-1; i >= distance; --i) {
                                if (cursor.equals(kmers.get(i)) &&
                                        partner.equals(kmers.get(i-distance))) {
                                    return false;
                                }
                            }
                        }

                        if (assembledKmersBloomFilter.lookup(cursor.getHash())) {
                            boolean assembled = greedyExtendLeft(graph, cursor, lookahead, lookahead, assembledKmersBloomFilter) != null;

                            if (assembled) {
    //                            int numNotAssembled = 0;

                                for (Kmer kmer : extension) {
                                    if (!assembledKmersBloomFilter.lookup(kmer.getHash())) {
    //                                    if (distanceInversePI < ++numNotAssembled) {
                                            assembled = false;
                                            break;
    //                                    }
                                    }
                                }

                                if (assembled) {
                                    for (int i=partnerIndex; i<numKmers; ++i) {
                                        if (!assembledKmersBloomFilter.lookup(kmers.get(i).getHash())) {
    //                                        if (distanceInversePI < ++numNotAssembled) {
                                                assembled = false;
                                                break;
    //                                        }
                                        }
                                    }
                                }
                            }

                            if (assembled) {
                                // region already assembled, do not extend further
                                return false;
                            }
                        }

                        kmers.addAll(extension);
                        kmers.add(cursor);

                        usedKmers.addAll(extension);
                        usedKmers.add(cursor);

                        return true;
                    }
                    else if (depth < maxDepth &&
                            (depth == 0 || !extensionKmers.contains(cursor))) {

                        if (cursor.hasAtLeastXSuccessors(k, numHash, graph, 2)) {
                            // only consider kmers that may be visited from an alternative branch upstream

                            ArrayDeque<Integer> visitedDepths = visitedKmers.get(cursor);
                            if (visitedDepths == null) {
                                visitedDepths = new ArrayDeque<>();
                                visitedDepths.add(depth);

                                visitedKmers.put(cursor, visitedDepths);

                                branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursor);
                                ++depth;
                                ++partnerIndex;
                            }
                            else {
                                boolean visited = false;
                                for (Integer d : visitedDepths) {
                                    if (Math.abs(d - depth) <= maxIndelSize) {
                                        visited = true;
                                        break;
                                    }
                                }

                                if (!visited) {
                                    visitedDepths.add(depth);

                                    branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                                    extension.add(cursor);
                                    extensionKmers.add(cursor);
                                    ++depth;
                                    ++partnerIndex;
                                }
                            }
                        }
                        else {
                            branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursor);
                            ++depth;
                            ++partnerIndex;
                        }
                    }
                }
            }
        }
        
        return false;
    }
    
    public static void extendWithPairedKmersBFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            float minKmerCov) {
        
        HashSet<Kmer> usedKmers = new HashSet<>(kmers);
        
        // naive extend LEFT
        Collections.reverse(kmers);
                
        // extend with paired kmers LEFT
        
        extendLeftWithPairedKmersBFS(kmers, 
                                    graph,
                                    maxTipLength,
                                    maxIndelSize,
                                    percentIdentity,
                                    minNumPairs,
                                    usedKmers,
                                    minKmerCov);
        
        Collections.reverse(kmers);
        
        // extend with paired kmers RIGHT
        
        extendRightWithPairedKmersBFS(kmers, 
                                    graph,
                                    maxTipLength,
                                    maxIndelSize,
                                    percentIdentity,
                                    minNumPairs,
                                    usedKmers,
                                    minKmerCov);
        
    }
    
    public static void extendWithPairedKmers(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            BloomFilter assembledKmersBloomFilter,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            float maxCovGradient,
                                            float minKmerCov) {
                
        HashSet<Kmer> usedKmers = new HashSet<>(kmers);
        
        // extend with paired kmers LEFT
        Collections.reverse(kmers);
                
        boolean extendable = extendLeftWithPairedKmersOnly(kmers, graph, usedKmers, assembledKmersBloomFilter);
        
        while (extendable) {
            extendable = extendLeftWithPairedKmersBFS(kmers, 
                                        graph, 
                                        lookahead, 
                                        maxTipLength,
                                        maxIndelSize,
                                        percentIdentity,
                                        minNumPairs,
                                        assembledKmersBloomFilter,
                                        usedKmers,
                                        maxCovGradient,
                                        minKmerCov);

//            System.out.println(graph.assembleReverseOrder(new ArrayDeque<>(kmers)));
            
            if (extendable) {
                extendable = extendLeftWithPairedKmersDFS(kmers, 
                                        graph, 
                                        lookahead, 
                                        maxTipLength,
                                        maxIndelSize,
                                        percentIdentity,
                                        minNumPairs,
                                        assembledKmersBloomFilter,
                                        usedKmers,
                                        maxCovGradient);
                
//                System.out.println(graph.assembleReverseOrder(new ArrayDeque<>(kmers)));
                
                if (extendable) {
                    extendable = extendLeftWithPairedKmersOnly(kmers, graph, usedKmers, assembledKmersBloomFilter);                    
                    
//                    System.out.println(graph.assembleReverseOrder(new ArrayDeque<>(kmers)));
                }
            }
        }
        
        Collections.reverse(kmers);
                
        // extend with paired kmers RIGHT
        
        extendable = extendRightWithPairedKmersOnly(kmers, graph, usedKmers, assembledKmersBloomFilter);
        
        while (extendable) {
            extendable = extendRightWithPairedKmersBFS(kmers, 
                                        graph, 
                                        lookahead, 
                                        maxTipLength,
                                        maxIndelSize,
                                        percentIdentity,
                                        minNumPairs,
                                        assembledKmersBloomFilter,
                                        usedKmers,
                                        maxCovGradient,
                                        minKmerCov);
            
            if (extendable) {
                extendable = extendRightWithPairedKmersDFS(kmers, 
                                        graph, 
                                        lookahead, 
                                        maxTipLength,
                                        maxIndelSize,
                                        percentIdentity,
                                        minNumPairs,
                                        assembledKmersBloomFilter,
                                        usedKmers,
                                        maxCovGradient);
                
                if (extendable) {
                    extendable = extendRightWithPairedKmersOnly(kmers, graph, usedKmers, assembledKmersBloomFilter);
                }
            }
        }
    }

    public static void extendWithPairedKmersDFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            BloomFilter assembledKmersBloomFilter,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            float maxCovGradient) {
        
        HashSet<Kmer> usedKmers = new HashSet<>(kmers);
        
        // extend with paired kmers LEFT
        Collections.reverse(kmers);
        
        boolean extendable = extendLeftWithPairedKmersOnly(kmers, graph, usedKmers, assembledKmersBloomFilter);
        
        while (extendable) {
            extendable = extendLeftWithPairedKmersDFS(kmers, 
                                    graph, 
                                    lookahead, 
                                    maxTipLength,
                                    maxIndelSize,
                                    percentIdentity,
                                    minNumPairs,
                                    assembledKmersBloomFilter,
                                    usedKmers,
                                    maxCovGradient);
            
            if (extendable) {
                extendable = extendLeftWithPairedKmersOnly(kmers, graph, usedKmers, assembledKmersBloomFilter);
            }
        }
        
        Collections.reverse(kmers);
        
        // extend with paired kmers RIGHT
        
        extendable = extendRightWithPairedKmersOnly(kmers, graph, usedKmers, assembledKmersBloomFilter);
        
        while (extendable) {
            extendable = extendRightWithPairedKmersDFS(kmers, 
                                    graph, 
                                    lookahead, 
                                    maxTipLength,
                                    maxIndelSize,
                                    percentIdentity,
                                    minNumPairs,
                                    assembledKmersBloomFilter,
                                    usedKmers,
                                    maxCovGradient);
            
            if (extendable) {
                extendable = extendRightWithPairedKmersOnly(kmers, graph, usedKmers, assembledKmersBloomFilter);
            }
        }
    }
    
/*
    public static void extendWithPairedKmers(ArrayList<Kmer2> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength, 
                                            boolean greedy, 
                                            BloomFilter assembledKmersBloomFilter,
                                            int maxIndelSize,
                                            float percentIdentity) {
        
        final int distance = graph.getPairedKmerDistance();
                
        // kmer pairs used in the extension of this transcript
        final HashSet<String> usedPairs = new HashSet<>();
        
        final HashSet<String> usedKmers = new HashSet<>();
        for (Kmer kmer : kmers) {
            usedKmers.add(bytesToString(kmer.bytes, k));
        }
        
//        // naive extend left
//        ArrayList<Kmer2> naiveExtension = naiveExtendLeft(kmers.get(0), graph, maxTipLength, usedKmers, true);
//        if (!naiveExtension.isEmpty()) {
//            kmers.addAll(0, naiveExtension);
//        }
//        
//        // naive extend right
//        naiveExtension = naiveExtendRight(kmers.get(kmers.size()-1), graph, maxTipLength, usedKmers);
//        if (!naiveExtension.isEmpty()) {
//            kmers.addAll(naiveExtension);
//        }
        
        // extend right with paired kmers
  
        Kmer2 kmer, n, n2, partner, partner2;
        int partnerIndex;
        String mergedSeq;
        
        ArrayDeque<Kmer2> extension = new ArrayDeque<>();
        
        ArrayDeque<LinkedList<Kmer2>> branchesStack = new ArrayDeque<>();
        LinkedList<Kmer2> neighbors = getSuccessorsRanked(kmers.get(kmers.size()-1), graph, lookahead);
        branchesStack.add(neighbors);
        boolean stop = false;
        int depth, numKmers;
        
        HashSet<String> visitedKmers = new HashSet<>();

        int maxRightDepth = maxRightPartnerSearchDepth(kmers, graph, distance, assembledKmersBloomFilter);
        
        while (!branchesStack.isEmpty() && !stop && maxRightDepth > 0) {
            neighbors = branchesStack.getLast();
            depth = extension.size();
            
            if (neighbors.isEmpty()) {
                branchesStack.removeLast();
//                    extension.removeLast().successors = null; // prune the cached successors
                extension.pollLast();
            }
            else if (depth >= maxRightDepth) {
                numKmers = kmers.size();
                partnerIndex = numKmers - distance + depth;
                partner = kmers.get(partnerIndex);

                boolean found = false;
                
                if (graph.lookupLeftKmer(partner.hashVals)) {
                    partner2 = kmers.get(partnerIndex-1);
                    
                    if (graph.lookupLeftKmer(partner2.hashVals)) {

                        n2 = extension.peekLast();
                        if (n2 == null) {
                            n2 = kmers.get(numKmers-1);
                        }

                        if (graph.lookupRightKmer(n2.hashVals) &&
                                graph.lookupKmerPairing(partner2.hashVals, n2.hashVals)) {

                            ListIterator<Kmer2> nItr = neighbors.listIterator();
                            while (nItr.hasNext()) {
                                n = nItr.next();
                                nItr.remove();
                                
                                if (graph.lookupRightKmer(n.hashVals) && 
                                        graph.lookupKmerPairing(partner.hashVals, n.hashVals)) {

                                    branchesStack.clear();
                                    visitedKmers.clear();
                                    
                                    if (!greedy && 
                                            assembledKmersBloomFilter.lookup(n.hashVals) &&
                                            assembledKmersBloomFilter.lookup(partner.hashVals) &&
                                            assembledKmersBloomFilter.lookup(partner2.hashVals) &&
                                            assembledKmersBloomFilter.lookup(n2.hashVals) ) {
                                        stop = true;
                                        break;
                                    }
                                    
                                    kmers.addAll(extension);
                                    kmers.add(n);
                                    
                                    Iterator<Kmer2> kmerItr = extension.iterator();
                                    while (kmerItr.hasNext()) {
                                        kmer = kmerItr.next();
        //                                kmer.successors = null; // prune the cached successors
                                        usedKmers.add(kmer.seq);
                                        kmerItr.remove();
                                    }
                                    usedKmers.add(n.seq);

                                    mergedSeq = partner.seq + n.seq;

                                    if (usedPairs.contains(mergedSeq)) {
                                        stop = true;
                                        break;
                                    }
                                    else {
                                        usedPairs.add(mergedSeq);
                                    }

                                    ArrayDeque<Kmer2> extendRight = extendRight(n, graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
                                    if (extendRight.isEmpty()) {
                                        branchesStack.add(getSuccessorsRanked(n, graph, lookahead));
                                    }
                                    else {
                                        kmers.addAll(extendRight);
                                        for (Kmer e : extendRight) {
                                            usedKmers.add(e.seq);
                                        }
                                        
                                        branchesStack.add(getSuccessorsRanked(kmers.get(kmers.size()-1), graph, lookahead));
                                    }
                                    
                                    maxRightDepth = maxRightPartnerSearchDepth(kmers, graph, distance, assembledKmersBloomFilter);
        //                            n.successors = null; // prune the cached successors

                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (stop) {
                    break;
                }
                else if (!found) {
                    neighbors.clear();
                }
            }
            else {
                
                n = neighbors.removeFirst();
                numKmers = kmers.size();
                
                if (numKmers + depth >= distance + 1) {
                    partnerIndex = numKmers - distance + depth;
                    partner = kmers.get(partnerIndex);
                    mergedSeq = partner.seq + n.seq;
                    
                    if (!visitedKmers.contains(mergedSeq)) {
                        partner2 = kmers.get(partnerIndex-1);
                        
                        n2 = extension.peekLast();                        
                        if (n2 == null) {
                            n2 = kmers.get(numKmers-1);
                        }

                        if (graph.lookupKmerPair(partner.hashVals, n.hashVals) &&
                                graph.lookupKmerPair(partner2.hashVals, n2.hashVals)) {
                            
                            branchesStack.clear();
                            visitedKmers.clear();
                                                        
                            if (!greedy && 
                                    assembledKmersBloomFilter.lookup(n.hashVals) &&
                                    assembledKmersBloomFilter.lookup(partner.hashVals) &&
                                    assembledKmersBloomFilter.lookup(partner2.hashVals) &&
                                    assembledKmersBloomFilter.lookup(n2.hashVals)) {
                                stop = true;
                                break;
                            }
                            
                            kmers.addAll(extension);
                            kmers.add(n);
                            
                            Iterator<Kmer2> kmerItr = extension.iterator();
                            while (kmerItr.hasNext()) {
                                kmer = kmerItr.next();
//                                kmer.successors = null; // prune the cached successors
                                usedKmers.add(kmer.seq);
                                kmerItr.remove();
                            }
                            usedKmers.add(n.seq);
                            
                            if (usedPairs.contains(mergedSeq)) {
                                stop = true;
                                break;
                            }
                            else {
                                usedPairs.add(mergedSeq);
                            }
                            
                            ArrayDeque<Kmer2> extendRight = extendRight(n, graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
                            if (extendRight.isEmpty()) {
                                branchesStack.add(getSuccessorsRanked(n, graph, lookahead));
                            }
                            else {
                                kmers.addAll(extendRight);
                                for (Kmer e : extendRight) {
                                    usedKmers.add(e.seq);
                                }
                                
                                branchesStack.add(getSuccessorsRanked(kmers.get(kmers.size()-1), graph, lookahead));
                            }
                            
                            maxRightDepth = maxRightPartnerSearchDepth(kmers, graph, distance, assembledKmersBloomFilter);
                            
//                            n.successors = null; // prune the cached successors
                        }
                        else {
                            extension.add(n);
                            branchesStack.add(getSuccessorsRanked(n, graph, lookahead));
                            visitedKmers.add(mergedSeq);
                        }
                    }
                }
                else {
                    extension.add(n);
                    branchesStack.add(getSuccessorsRanked(n, graph, lookahead));
                }
            }
        }
        
        if (!stop) {
            // naive extend right for the final time
            ArrayDeque<Kmer2> naiveExtension = extendRight(kmers.get(kmers.size()-1), graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
            if (!naiveExtension.isEmpty()) {
                kmers.addAll(naiveExtension);
            }
        }
        
        // extend left with paired kmers
        
        Collections.reverse(kmers);
        
        stop = false;
        visitedKmers.clear();
        extension.clear();
        branchesStack.clear();
        neighbors = getPredecessorsRanked(kmers.get(kmers.size()-1), graph, lookahead);
        branchesStack.add(neighbors);
        int maxLeftDepth = maxLeftPartnerSearchDepth(kmers, graph, distance, assembledKmersBloomFilter);
        
        while (!branchesStack.isEmpty() && !stop && maxLeftDepth > 0) {
            neighbors = branchesStack.getLast();
            depth = extension.size();
            
            if (neighbors.isEmpty()) {
                branchesStack.removeLast();
//                    extension.removeLast().predecessors = null; // prune the cached predecessors
                extension.pollLast();
            }
            else if (depth >= maxLeftDepth) {
                numKmers = kmers.size();
                partnerIndex = numKmers - distance + depth;
                partner = kmers.get(partnerIndex);
                
                boolean found = false;
                
                if (graph.lookupRightKmer(partner.hashVals)) {
                    partner2 = kmers.get(partnerIndex-1);
                                    
                    if (graph.lookupRightKmer(partner2.hashVals)) {
                    
                        n2 = extension.peekLast();
                        
                        if (n2 == null) {
                            n2 = kmers.get(numKmers-1);
                        }

                        if (graph.lookupLeftKmer(n2.hashVals) && 
                                graph.lookupKmerPairing(n2.hashVals, partner2.hashVals)) {

                            ListIterator<Kmer2> nItr = neighbors.listIterator();
                            while (nItr.hasNext()) {
                                n = nItr.next();
                                nItr.remove();

                                if (graph.lookupLeftKmer(n.hashVals) && 
                                        graph.lookupKmerPairing(n.hashVals, partner.hashVals)) {

                                    branchesStack.clear();
                                    visitedKmers.clear();
                                    
                                    if (!greedy && 
                                            assembledKmersBloomFilter.lookup(n.hashVals) &&
                                            assembledKmersBloomFilter.lookup(partner.hashVals) &&
                                            assembledKmersBloomFilter.lookup(partner2.hashVals) &&
                                            assembledKmersBloomFilter.lookup(n2.hashVals)) {
                                        stop = true;
                                        break;
                                    }
                                    
                                    kmers.addAll(extension);
                                    kmers.add(n);
                                    
                                    Iterator<Kmer2> kmerItr = extension.iterator();
                                    while (kmerItr.hasNext()) {
                                        kmer = kmerItr.next();
        //                                kmer.predecessors = null; // prune the cached predecessors
                                        usedKmers.add(kmer.seq);
                                        kmerItr.remove();
                                    }
                                    usedKmers.add(n.seq);

                                    mergedSeq = n.seq + partner.seq;

                                    if (usedPairs.contains(mergedSeq)) {
                                        stop = true;
                                        break;
                                    }
                                    else {
                                        usedPairs.add(mergedSeq);
                                    }

                                    ArrayDeque<Kmer2> extendLeft = extendLeft(n, graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
                                    if (extendLeft.isEmpty()) {
                                        branchesStack.add(getPredecessorsRanked(n, graph, lookahead));
                                    }
                                    else {
                                        kmers.addAll(extendLeft);
                                        for (Kmer e : extendLeft) {
                                            usedKmers.add(e.seq);
                                        }
                                        
                                        branchesStack.add(getPredecessorsRanked(kmers.get(kmers.size()-1), graph, lookahead));
                                    }
                                    
                                    maxLeftDepth = maxLeftPartnerSearchDepth(kmers, graph, distance, assembledKmersBloomFilter);

        //                            n.predecessors = null; // prune the cached predecessors

                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                
                if (stop) {
                    break;
                }
                else if (!found) {
                    neighbors.clear();
                }
            }
            else {
                
                n = neighbors.removeFirst();
                numKmers = kmers.size();
                
                if (numKmers + depth >= distance + 1) {
                    partnerIndex = numKmers - distance + depth;
                    partner = kmers.get(partnerIndex);
                    mergedSeq = n.seq + partner.seq;

                    if (!visitedKmers.contains(mergedSeq)) {
                        partner2 = kmers.get(partnerIndex-1);

                        n2 = extension.peekLast();
                        if (n2 == null) {
                            n2 = kmers.get(numKmers-1);
                        }
                        
                        if (graph.lookupKmerPair(n.hashVals, partner.hashVals) &&
                                graph.lookupKmerPair(n2.hashVals, partner2.hashVals)) {
                            
                            branchesStack.clear();
                            visitedKmers.clear();
                            
                            if (!greedy && 
                                    assembledKmersBloomFilter.lookup(n.hashVals) &&
                                    assembledKmersBloomFilter.lookup(partner.hashVals) &&
                                    assembledKmersBloomFilter.lookup(partner2.hashVals) &&
                                    assembledKmersBloomFilter.lookup(n2.hashVals)) {
                                stop = true;
                                break;
                            }
                            
                            kmers.addAll(extension);
                            kmers.add(n);
                            
                            Iterator<Kmer2> kmerItr = extension.iterator();
                            while (kmerItr.hasNext()) {
                                kmer = kmerItr.next();
//                                kmer.predecessors = null; // prune the cached predecessors
                                usedKmers.add(kmer.seq);
                                kmerItr.remove();
                            }
                            usedKmers.add(n.seq);

                            if (usedPairs.contains(mergedSeq)) {
                                stop = true;
                            }
                            else {
                                usedPairs.add(mergedSeq);
                            }

                            ArrayDeque<Kmer2> extendLeft = extendLeft(n, graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
                            if (extendLeft.isEmpty()) {
                                branchesStack.add(getPredecessorsRanked(n, graph, lookahead));
                            }
                            else {
                                kmers.addAll(extendLeft);
                                for (Kmer e : extendLeft) {
                                    usedKmers.add(e.seq);
                                }
                                
                                branchesStack.add(getPredecessorsRanked(kmers.get(kmers.size()-1), graph, lookahead));
                            }
                            
                            maxLeftDepth = maxLeftPartnerSearchDepth(kmers, graph, distance, assembledKmersBloomFilter);
//                            n.predecessors = null; // prune the cached predecessors
                        }
                        else {
                            extension.add(n);
                            branchesStack.add(getPredecessorsRanked(n, graph, lookahead));
                            visitedKmers.add(mergedSeq);
                        }
                    }
                }
                else {
                    extension.add(n);
                    branchesStack.add(getPredecessorsRanked(n, graph, lookahead));
                }
            }
        }
        
        if (!stop) {
            // naive extend left for the final time
            // note: list is still reversed here
            ArrayDeque<Kmer2> naiveExtension = extendLeft(kmers.get(kmers.size()-1), graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
            if (!naiveExtension.isEmpty()) {
                kmers.addAll(naiveExtension);
            }
        }
        
        Collections.reverse(kmers);
        
        //return assemble(kmers, k);
    }
*/
    
//    private static boolean hasFragmentDepthRight(String source, BloomFilterDeBruijnGraph graph, int depth) {
//        LinkedList<LinkedList> frontier = new LinkedList<>();
//        LinkedList<String> alts = new LinkedList<>();
//        for (String s : graph.getSuccessors(source)) {
//            if (graph.lookupFragmentKmer(s)) {
//                alts.add(s);
//            }
//        }
//        frontier.add(alts);
//        
//        while (!frontier.isEmpty()) {
//            alts = frontier.peekLast();
//            if (alts.isEmpty()) {
//                frontier.removeLast();
//            }
//            else {
//                String a = alts.pop();
//                alts = new LinkedList<>();
//                for (String s : graph.getSuccessors(a)) {
//                    if (graph.lookupFragmentKmer(s)) {
//                        alts.add(s);
//                    }
//                }
//                frontier.add(alts);
//            }
//
//            if (frontier.size() >= depth) {
//                return true;
//            }
//        }
//        
//        return false;
//    }
//    
//    private static boolean hasFragmentDepthLeft(String source, BloomFilterDeBruijnGraph graph, int depth) {
//        LinkedList<LinkedList> frontier = new LinkedList<>();
//        LinkedList<String> alts = new LinkedList<>();
//        for (String s : graph.getPredecessors(source)) {
//            if (graph.lookupFragmentKmer(s)) {
//                alts.add(s);
//            }
//        }
//        frontier.add(alts);
//        
//        while (!frontier.isEmpty()) {
//            alts = frontier.peekLast();
//            if (alts.isEmpty()) {
//                frontier.removeLast();
//            }
//            else {
//                String a = alts.pop();
//                alts = new LinkedList<>();
//                for (String s : graph.getPredecessors(a)) {
//                    if (graph.lookupFragmentKmer(s)) {
//                        alts.add(s);
//                    }
//                }
//                frontier.add(alts);
//            }
//
//            if (frontier.size() >= depth) {
//                return true;
//            }
//        }
//        
//        return false;
//    }
        
    public static boolean hasDepthRight(Kmer source, BloomFilterDeBruijnGraph graph, int depth) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<ArrayDeque> frontier = new ArrayDeque<>();
        ArrayDeque<Kmer> alts = source.getSuccessors(k, numHash, graph);
        frontier.add(alts);
        
        while (!frontier.isEmpty()) {
            alts = frontier.peekLast();
            if (alts.isEmpty()) {
                frontier.removeLast();
            }
            else {
                frontier.add(alts.pop().getSuccessors(k, numHash, graph));
            }

            if (frontier.size() >= depth) {
                return true;
            }
        }
        
        return false;
    }
    
    public static boolean hasDepthLeft(Kmer source, BloomFilterDeBruijnGraph graph, int depth) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<ArrayDeque> frontier = new ArrayDeque<>();
        ArrayDeque<Kmer> alts = source.getPredecessors(k, numHash, graph);
        frontier.add(alts);
        
        while (!frontier.isEmpty()) {
            alts = frontier.peekLast();
            if (alts.isEmpty()) {
                frontier.removeLast();
            }
            else {
                frontier.add(alts.pop().getPredecessors(k, numHash, graph));
            }

            if (frontier.size() >= depth) {
                return true;
            }
        }
        
        return false;
    }
    
    public static boolean hasDepthRight(Kmer source, BloomFilterDeBruijnGraph graph, int depth, BloomFilter bf) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<ArrayDeque> frontier = new ArrayDeque<>();
        ArrayDeque<Kmer> alts = source.getSuccessors(k, numHash, graph, bf);
        frontier.add(alts);
        
        while (!frontier.isEmpty()) {
            alts = frontier.peekLast();
            if (alts.isEmpty()) {
                frontier.removeLast();
            }
            else {
                frontier.add(alts.pop().getSuccessors(k, numHash, graph, bf));
            }

            if (frontier.size() >= depth) {
                return true;
            }
        }
        
        return false;
    }
    
    public static boolean hasDepthLeft(Kmer source, BloomFilterDeBruijnGraph graph, int depth, BloomFilter bf) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<ArrayDeque> frontier = new ArrayDeque<>();
        ArrayDeque<Kmer> alts = source.getPredecessors(k, numHash, graph, bf);
        frontier.add(alts);
        
        while (!frontier.isEmpty()) {
            alts = frontier.peekLast();
            if (alts.isEmpty()) {
                frontier.removeLast();
            }
            else {
                frontier.add(alts.pop().getPredecessors(k, numHash, graph, bf));
            }

            if (frontier.size() >= depth) {
                return true;
            }
        }
        
        return false;
    }
    
    public static ArrayDeque<Kmer> naiveExtendRight(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<Kmer> terminators, float minKmerCov) {        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        HashSet<Kmer> usedKmers = new HashSet<>();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        kmer.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : best.getLeftVariants(k, numHash, graph)) {
                if (s.hasDepthLeft(k, numHash, graph, maxTipLength)) {
                    return result;
                }
            }
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
            }
            else {
                best = null;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    if (n.hasDepthRight(k, numHash, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            return result;
                        }
                    }
                }
            }
            
            if (best == null) {
                break;
            }
                        
            if (terminators.contains(best) || usedKmers.contains(best)) {
                break;
            }
            
            result.add(best);
            usedKmers.add(best);
            
            best.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> naiveExtendRight(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, int bound, float minKmerCov) {        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        int extensionLength = 0;
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        kmer.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : best.getLeftVariants(k, numHash, graph)) {
                if (s.hasDepthLeft(k, numHash, graph, maxTipLength)) {
                    return result;
                }
            }
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
            }
            else {
                best = null;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    if (n.hasDepthRight(k, numHash, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            return result;
                        }
                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            result.add(best);
            
            if (++extensionLength > bound) {
                break;
            }
                        
            best.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> naiveExtendRightNoBackChecks(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, int bound, float minKmerCov) {        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        int extensionLength = 0;
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        kmer.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        while (!neighbors.isEmpty()) {
            Kmer best;
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
            }
            else {
                best = null;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    if (n.hasDepthRight(k, numHash, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            return result;
                        }
                    }
                }
            }
            
            if (best == null || best.equals(kmer) || (!result.isEmpty() && best.equals(result.peekLast()))) {
                break;
            }
            
            result.add(best);
            
            if (++extensionLength > bound) {
                break;
            }
                        
            best.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayList<Kmer> naiveExtend(ArrayList<Kmer> kmers, BloomFilterDeBruijnGraph graph, int maxTipLength, float minKmerCov) {
        HashSet<Kmer> usedKmers = new HashSet<>(kmers);
        
        ArrayDeque<Kmer> leftExtension = naiveExtendLeft(kmers.get(0), graph, maxTipLength, usedKmers, minKmerCov);
        usedKmers.addAll(leftExtension);
        
        ArrayDeque<Kmer> rightExtension = naiveExtendRight(kmers.get(kmers.size()-1), graph, maxTipLength, usedKmers, minKmerCov);
        
        ArrayList<Kmer> result = new ArrayList<>(leftExtension.size() + kmers.size() + rightExtension.size());
        
        if (!leftExtension.isEmpty()) {
            Iterator<Kmer> itr = leftExtension.descendingIterator();
            while (itr.hasNext()) {
                result.add(itr.next());
            }
        }
        
        result.addAll(kmers);
        
        result.addAll(rightExtension);
        
        return result;
    }
    
    public static ArrayDeque<Kmer> naiveExtendLeft(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<Kmer> terminators, float minKmerCov) {        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        HashSet<Kmer> usedKmers = new HashSet<>();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        kmer.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : best.getRightVariants(k, numHash, graph)) {
                if (s.hasDepthRight(k, numHash, graph, maxTipLength)) {
                    return result;
                }
            }
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
            }
            else {
                best = null;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    if (n.hasDepthLeft(k, numHash, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            return result;
                        }
                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            if (terminators.contains(best) || usedKmers.contains(best)) {
                break;
            }
            
            result.addLast(best);
            usedKmers.add(best);
            
            best.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> naiveExtendLeft(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, int bound, float minKmerCov) {        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        int extensionLength = 0;
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        kmer.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : best.getRightVariants(k, numHash, graph)) {
                if (s.hasDepthRight(k, numHash, graph, maxTipLength)) {
                    return result;
                }
            }
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
            }
            else {
                best = null;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    if (n.hasDepthLeft(k, numHash, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            return result;
                        }
                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            result.addLast(best);
            
            if (++extensionLength > bound) {
                break;
            }
            
            best.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> naiveExtendLeftNoBackChecks(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, int bound, float minKmerCov) {        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        int extensionLength = 0;
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        kmer.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        while (!neighbors.isEmpty()) {
            Kmer best;
                    
            if (neighbors.size() == 1) {
                best = neighbors.pop();
            }
            else {
                best = null;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    if (n.hasDepthLeft(k, numHash, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            return result;
                        }
                    }
                }
            }
            
            if (best == null || best.equals(kmer) || (!result.isEmpty() && best.equals(result.peekLast()))) {
                break;
            }
            
            result.addLast(best);
            
            if (++extensionLength > bound) {
                break;
            }
            
            best.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> extendRight(Kmer source,
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLength,
                                            int bound, 
                                            int maxIndelSize, 
                                            float percentIdentity,
                                            float minKmerCov) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        source.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        Kmer best;
        
        while (!neighbors.isEmpty() && result.size() <= bound) {
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
                result.add(best);
                
                ArrayDeque<Kmer> b = naiveExtendRight(best, graph, maxTipLength, bound-result.size(), minKmerCov);
                
                if (b.isEmpty()) {
                    break;
                }
                else {
                    result.addAll(b);
                    best = result.peekLast();
                }
            }
            else {
                best = null;
                
                ArrayDeque<ArrayDeque<Kmer>> branches = new ArrayDeque<>(4);
                
                float maxCov = -1;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();

                    if (n.hasDepthRight(k, numHash, graph, maxTipLength)) {
                        ArrayDeque<Kmer> b = naiveExtendRight(n, graph, maxTipLength, bound-result.size(), minKmerCov);
                        
                        if (b.size() < maxTipLength) {
                            // indicates too many branches; can't resolve
                            return result;
                        }
                        
                        b.addFirst(n);
                        
                        float c = getMedianKmerCoverage(b);
                        if (c > maxCov) {
                            maxCov = c;
                            branches.addFirst(b);
                        }
                        else {
                            branches.addLast(b);
                        }
                    }
                }
                
                if (branches.isEmpty()) {
                    // fuzzy end
                    return result;
                }
                else if (branches.size() == 1) {
                    // one good branch
                    
                    ArrayDeque<Kmer> b = branches.pop();
                    result.addAll(b);
                    best = result.peekLast();
                }
                else {
                    // multiple good branches
                    
                    // check whether the branches form a bubble
                    ArrayDeque<Kmer> bestBranch = branches.pollFirst();
                    String bestBranchSeq = graph.assemble(bestBranch);
                    String suffix = graph.getSuffix(bestBranch.peekLast().toString());
                    int bestBranchLength = bestBranch.size();
                    
                    for (ArrayDeque<Kmer> b : branches) {
                        int len = b.size();
                        
                        if (len >= bestBranchLength - maxIndelSize && len <= bestBranchLength + maxIndelSize) {
                            // length is within range
                            
                            if (!suffix.equals(graph.getSuffix(b.peekLast().toString())) || 
                                    getPercentIdentity(graph.assemble(b), bestBranchSeq) < percentIdentity) {
                                return result;
                            }
                        }
                        else {
                            if (len < bestBranchLength - maxIndelSize) {
                                // compare percent identity
                                if (getPercentIdentity(graph.assemble(b), bestBranchSeq.substring(0, len+k-1)) < percentIdentity) {
                                    return result;
                                }
                            }
                            else {
                                return result;
                            }
                        }
                    }
                    
                    result.addAll(bestBranch);
                    best = bestBranch.peekLast();
                    best.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
                    if (neighbors.size() == 1) {
                        // bubble branches converge at this kmer
                        best = neighbors.pop();
                        result.add(best);
                    }
                    else {
                        neighbors.clear();
                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            best.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> extendLeft(Kmer source,
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLength, 
                                            int bound, 
                                            int maxIndelSize, 
                                            float percentIdentity,
                                            float minKmerCov) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        source.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        Kmer best;
        
        while (!neighbors.isEmpty()) {
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
                result.add(best);
                
                ArrayDeque<Kmer> b = naiveExtendLeft(best, graph, maxTipLength, bound-result.size(), minKmerCov);
                
                if (b.isEmpty()) {
                    break;
                }
                else {
                    result.addAll(b);
                    best = result.peekLast();
                }
            }
            else {
                best = null;
                
                ArrayDeque<ArrayDeque<Kmer>> branches = new ArrayDeque<>(4);
                
                float maxCov = -1;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    
//                    if (hasDepthLeft(n, graph, maxTipLength)) {
                    if (n.hasDepthLeft(k, numHash, graph, maxTipLength)) {
                        ArrayDeque<Kmer> b = naiveExtendLeft(n, graph, maxTipLength, bound-result.size(), minKmerCov);
                        if (b.size() < maxTipLength) {
                            // indicates too many branches; can't resolve
                            return result;
                        }
                        
                        b.addFirst(n);
                        
                        float c = getMedianKmerCoverage(b);
                        if (c > maxCov) {
                            maxCov = c;
                            branches.addFirst(b);
                        }
                        else {
                            branches.addLast(b);
                        }
                    }
                }
                
                if (branches.isEmpty()) {
                    // fuzzy end
                    return result;
                }
                else if (branches.size() == 1) {
                    // one good branch
                    ArrayDeque<Kmer> b = branches.pop();
                    result.addAll(b);
                    best = result.peekLast();
                }
                else {
                    // multiple good branches
                    
                    // check whether the branches form a bubble
                    ArrayDeque<Kmer> bestBranch = branches.pollFirst();
                    String bestBranchSeq = graph.assembleReverseOrder(bestBranch);
                    String prefix = graph.getPrefix(bestBranch.peekLast().toString());
                    int bestBranchLength = bestBranch.size();
                    
                    for (ArrayDeque<Kmer> b : branches) {
                        int len = b.size();
                        
                        if (len >= bestBranchLength - maxIndelSize && len <= bestBranchLength + maxIndelSize) {
                            // length is within range
                            
                            if (!prefix.equals(graph.getPrefix(b.peekLast().toString())) || 
                                    getPercentIdentity(graph.assembleReverseOrder(b), bestBranchSeq) < percentIdentity) {
                                return result;
                            }
                        }
                        else {
                            if (len < bestBranchLength - maxIndelSize) {
                                // compare percent identity
                                if (getPercentIdentity(graph.assembleReverseOrder(b), bestBranchSeq.substring(bestBranchSeq.length() - (len+k-1))) < percentIdentity) {
                                    return result;
                                }
                            }
                            else {
                                return result;
                            }
                        }
                    }
                    
                    result.addAll(bestBranch);
                    best = bestBranch.peekLast();
                    best.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
                    if (neighbors.size() == 1) {
                        // bubble branches converge at this kmer
                        best = neighbors.pop();
                        result.add(best);
                    }
                    else {
                        neighbors.clear();
                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            best.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> extendRight(Kmer source,
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLength, 
                                            HashSet<Kmer> usedKmers, 
                                            int maxIndelSize, 
                                            float percentIdentity,
                                            float minKmerCov) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        source.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        Kmer best;
        
        while (!neighbors.isEmpty()) {
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
                
                if (usedKmers.contains(best)) {
                    return result;
                }
                else {
                    result.add(best);
                    usedKmers.add(best);
                }
                
                ArrayDeque<Kmer> b = naiveExtendRight(best, graph, maxTipLength, usedKmers, minKmerCov);
                
                if (b.isEmpty()) {
                    break;
                }
                else {
                    result.addAll(b);
                    usedKmers.addAll(b);
                    best = result.peekLast();
                }
            }
            else {
                best = null;
                
                ArrayDeque<ArrayDeque<Kmer>> branches = new ArrayDeque<>(4);
                
                float maxCov = -1;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    
//                    if (hasDepthRight(n, graph, maxTipLength)) {
                    if (n.hasDepthRight(k, numHash, graph, maxTipLength)) {
                        ArrayDeque<Kmer> b = naiveExtendRight(n, graph, maxTipLength, usedKmers, minKmerCov);
                        
                        if (b.size() < maxTipLength) {
                            // indicates too many branches; can't resolve
                            return result;
                        }
                        
                        b.addFirst(n);
                        
                        float c = getMedianKmerCoverage(b);
                        if (c > maxCov) {
                            maxCov = c;
                            branches.addFirst(b);
                        }
                        else {
                            branches.addLast(b);
                        }
                    }
                }
                
                if (branches.isEmpty()) {
                    // fuzzy end
                    return result;
                }
                else if (branches.size() == 1) {
                    // one good branch
                    
                    ArrayDeque<Kmer> b = branches.pop();
                    result.addAll(b);
                    usedKmers.addAll(b);
                    best = result.peekLast();
                }
                else {
                    // multiple good branches
                    
                    // check whether the branches form a bubble
                    ArrayDeque<Kmer> bestBranch = branches.pollFirst();
                    String bestBranchSeq = graph.assemble(bestBranch);
                    String suffix = graph.getSuffix(bestBranch.peekLast().toString());
                    int bestBranchLength = bestBranch.size();
                    
                    for (ArrayDeque<Kmer> b : branches) {
                        int len = b.size();
                        
                        if (len >= bestBranchLength - maxIndelSize && len <= bestBranchLength + maxIndelSize) {
                            // length is within range
                            
                            if (!suffix.equals(graph.getSuffix(b.peekLast().toString())) || 
                                    getPercentIdentity(graph.assemble(b), bestBranchSeq) < percentIdentity) {
                                return result;
                            }
                        }
                        else {
                            if (len < bestBranchLength - maxIndelSize) {
                                // compare percent identity
                                if (getPercentIdentity(graph.assemble(b), bestBranchSeq.substring(0, len+k-1)) < percentIdentity) {
                                    return result;
                                }
                            }
                            else {
                                return result;
                            }
                        }
                    }
                    
                    result.addAll(bestBranch);
                    usedKmers.addAll(bestBranch);
                    best = bestBranch.peekLast();
                    best.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
                    if (neighbors.size() == 1) {
                        // bubble branches converge at this kmer
                        best = neighbors.pop();
                        result.add(best);
                        usedKmers.add(best);
                    }
                    else {
                        neighbors.clear();
                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            best.getSuccessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> extendLeft(Kmer source,
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLength, 
                                            HashSet<Kmer> usedKmers, 
                                            int maxIndelSize, 
                                            float percentIdentity,
                                            float minKmerCov) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        source.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        Kmer best;
        
        while (!neighbors.isEmpty()) {
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
                
                if (usedKmers.contains(best)) {
                    return result;
                }
                
                result.add(best);
                
                ArrayDeque<Kmer> b = naiveExtendLeft(best, graph, maxTipLength, usedKmers, minKmerCov);
                
                if (b.isEmpty()) {
                    break;
                }
                else {
                    result.addAll(b);
                    usedKmers.addAll(b);
                    best = result.peekLast();
                }
            }
            else {
                best = null;
                
                ArrayDeque<ArrayDeque<Kmer>> branches = new ArrayDeque<>(4);
                
                float maxCov = -1;
                while (!neighbors.isEmpty()) {
                    Kmer n = neighbors.pop();
                    
//                    if (hasDepthLeft(n, graph, maxTipLength)) {
                    if (n.hasDepthLeft(k, numHash, graph, maxTipLength)) {
                        ArrayDeque<Kmer> b = naiveExtendLeft(n, graph, maxTipLength, usedKmers, minKmerCov);
                        if (b.size() < maxTipLength) {
                            // indicates too many branches; can't resolve
                            return result;
                        }
                        
                        b.addFirst(n);
                        
                        float c = getMedianKmerCoverage(b);
                        if (c > maxCov) {
                            maxCov = c;
                            branches.addFirst(b);
                        }
                        else {
                            branches.addLast(b);
                        }
                    }
                }
                
                if (branches.isEmpty()) {
                    // fuzzy end
                    return result;
                }
                else if (branches.size() == 1) {
                    // one good branch
                    ArrayDeque<Kmer> b = branches.pop();
                    result.addAll(b);
                    usedKmers.addAll(b);
                    best = result.peekLast();
                }
                else {
                    // multiple good branches
                    
                    // check whether the branches form a bubble
                    ArrayDeque<Kmer> bestBranch = branches.pollFirst();
                    String bestBranchSeq = graph.assembleReverseOrder(bestBranch);
                    String prefix = graph.getPrefix(bestBranch.peekLast().toString());
                    int bestBranchLength = bestBranch.size();
                    
                    for (ArrayDeque<Kmer> b : branches) {
                        int len = b.size();
                        
                        if (len >= bestBranchLength - maxIndelSize && len <= bestBranchLength + maxIndelSize) {
                            // length is within range
                            
                            if (!prefix.equals(graph.getPrefix(b.peekLast().toString())) || 
                                    getPercentIdentity(graph.assembleReverseOrder(b), bestBranchSeq) < percentIdentity) {
                                return result;
                            }
                        }
                        else {
                            if (len < bestBranchLength - maxIndelSize) {
                                // compare percent identity
                                if (getPercentIdentity(graph.assembleReverseOrder(b), bestBranchSeq.substring(bestBranchSeq.length() - (len+k-1))) < percentIdentity) {
                                    return result;
                                }
                            }
                            else {
                                return result;
                            }
                        }
                    }
                    
                    result.addAll(bestBranch);
                    usedKmers.addAll(bestBranch);
                    best = bestBranch.peekLast();
                    best.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
                    if (neighbors.size() == 1) {
                        // bubble branches converge at this kmer
                        best = neighbors.pop();
                        result.add(best);
                        usedKmers.add(best);
                    }
                    else {
                        neighbors.clear();
                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            best.getPredecessors(k, numHash, graph, neighbors, minKmerCov);
        }
        
        return result;
    }
    
    public static boolean isBranchFree(ArrayList<Kmer> seqKmers, BloomFilterDeBruijnGraph graph, int maxTipLength) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        for (Kmer kmer : seqKmers) {
            for (Kmer var : kmer.getRightVariants(k, numHash, graph)) {
//                if (hasDepthRight(var, graph, maxTipLength)) {
                if (var.hasDepthRight(k, numHash, graph, maxTipLength)) {
                    return false;
                }
            }
            
            for (Kmer var : kmer.getLeftVariants(k, numHash, graph)) {
//                if (hasDepthLeft(var, graph, maxTipLength)) {
                if (var.hasDepthLeft(k, numHash, graph, maxTipLength)) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    public static boolean isChimera(ArrayList<Kmer> seqKmers, BloomFilterDeBruijnGraph graph, BloomFilter assembledKmers, int lookahead) {
        int k = graph.getK();   
        int numKmers = seqKmers.size();
        int maxGap = 2*k;
        
        if (assembledKmers.lookup(seqKmers.get(0).getHash()) &&
                assembledKmers.lookup(seqKmers.get(numKmers-1).getHash())) {
            
            int i = 1;
            for (; i < numKmers-1; ++i) {
                if (!assembledKmers.lookup(seqKmers.get(i).getHash())) {
                    Kmer left = seqKmers.get(i-1);
                    
                    // check to see if this is a small gap
                    int t=i+1;
                    for (; t <numKmers-1; ++t) {
                        if (assembledKmers.lookup(seqKmers.get(t).getHash())) {
                            break;
                        }
                    }
                    
                    if (t < numKmers-1) {
                        Kmer right = seqKmers.get(t);
                        int d = t-i;
                        if (d <= maxGap && getMaxCoveragePath(graph, left, right, d, lookahead, assembledKmers) != null) {
                            i = t;
                            continue;
                        }
                    }
                    
                    break;
                }
            }
            
            if (i == numKmers -1) {
                return false;
            }
            
            --i;
            
            int j = numKmers-2;
            for (; j > i; --j) {
                if (!assembledKmers.lookup(seqKmers.get(j).getHash())) {
                    Kmer right = seqKmers.get(j+1);
                    
                    // check to see if this is a small gap
                    int t=j-1;
                    for (; t > i; --t) {
                        if (assembledKmers.lookup(seqKmers.get(t).getHash())) {
                            break;
                        }
                    }
                    
                    if (t > i) { 
                        Kmer left = seqKmers.get(t);
                        int d = j-t;
                        if (d <= maxGap && getMaxCoveragePath(graph, left, right, d, lookahead, assembledKmers) != null) {
                            j = t;
                            continue;
                        }
                    }
                    
                    break;
                }
            }
            
            ++j;
            
            if (j - i <= maxGap) { // use 2*k instead of k to account for mismatches near breakpoint
                Kmer source1 = seqKmers.get(i);
                Kmer source2 = seqKmers.get(j);
                
                HashSet<Kmer> kmers1 = new HashSet<>(greedyExtendRight(graph, source1, lookahead, 1000, assembledKmers));
                
                kmers1.retainAll(greedyExtendLeft(graph, source2, lookahead, 1000, assembledKmers));
                
                if (kmers1.isEmpty()) {
                    // 2 non-intersecting paths; most likely 2 transcripts from different genes
                    return true;
                }
                
                // otherwise, this is alternative splicing
            }
        }
        
        return false;
    }

    public static ArrayList<Kmer> trimReverseComplementArtifact(ArrayList<Kmer> seqKmers,
            BloomFilterDeBruijnGraph graph, boolean stranded,
            int maxEdgeClip, int maxIndelSize, float minPercentIdentity, float maxCovGradient) {
        
        int minMatchLen = 2*graph.getK();
        int numKmers = seqKmers.size();
        
        // scan from left to right
        int anchorStartIndex = -1;
        int anchorEndIndex = -1;
        int partnerStartIndex = -1;
        int partnerEndIndex = -1;
        int anchorStopIndex = Math.min(numKmers, maxEdgeClip);
        for (int i=0; i<anchorStopIndex; ++i) {
            Kmer seed = seqKmers.get(i);
            long rcHashVal = seed.getReverseComplementHash();
            for (int j=i+1; j<numKmers; ++j) {
                Kmer candidate = seqKmers.get(j);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    anchorStartIndex = i;
                    anchorEndIndex = i;
                    partnerStartIndex = j;
                    partnerEndIndex = j;
                    break;
                }
            }
            
            if (anchorStartIndex > 0) {
                break;
            }
        }
        
        if (anchorStartIndex > 0) {
            int midPointIndex = (anchorStartIndex + partnerStartIndex)/2;
            for (int i=anchorEndIndex+1; i<midPointIndex; ++i) {
                Kmer seed = seqKmers.get(i);
                long rcHashVal = seed.getReverseComplementHash();

                for (int j=partnerStartIndex-1; j>=midPointIndex; --j) {
                    Kmer candidate = seqKmers.get(j);
                    if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                        anchorEndIndex = i;
                        partnerStartIndex = j;
                        break;
                    }
                }
            }

            if (anchorStartIndex > 0 && anchorEndIndex-anchorStartIndex >= minMatchLen && partnerEndIndex-partnerStartIndex >= minMatchLen) {
                ++anchorEndIndex;
                ++partnerEndIndex;

                if (stranded) {
                    float anchorCov = getMedianKmerCoverage(seqKmers, anchorStartIndex, anchorEndIndex);
                    float middleCov = anchorEndIndex < partnerStartIndex ? getMedianKmerCoverage(seqKmers, anchorEndIndex, partnerStartIndex) : 0;
                    float partnerCov = getMedianKmerCoverage(seqKmers, partnerStartIndex, partnerEndIndex);

                    if (anchorCov < partnerCov) {
                        if (middleCov >= anchorCov && middleCov >= partnerCov * maxCovGradient) {
                            return new ArrayList<>(seqKmers.subList(anchorEndIndex, numKmers));
                        }
                        else {
                            return new ArrayList<>(seqKmers.subList(partnerStartIndex, numKmers));
                        }
                    }
                    else {
                        if (middleCov > partnerCov && middleCov >= anchorCov * maxCovGradient) {
                            return new ArrayList<>(seqKmers.subList(0, partnerStartIndex));
                        }
                        else {
                            return new ArrayList<>(seqKmers.subList(0, anchorEndIndex));
                        }
                    }
                }
                else {
                    return new ArrayList<>(seqKmers.subList(anchorEndIndex,partnerStartIndex));
                }
            }
        }
        
        // scan from right to left
        anchorStartIndex = -1;
        anchorEndIndex = -1;
        partnerStartIndex = -1;
        partnerEndIndex = -1;
        anchorStopIndex = Math.max(0, numKmers-maxEdgeClip);
        for (int i=numKmers-1; i>=anchorStopIndex; --i) {
            Kmer seed = seqKmers.get(i);
            long rcHashVal = seed.getReverseComplementHash();
            for (int j=i-1; j>=0; --j) {
                Kmer candidate = seqKmers.get(j);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    anchorStartIndex = i;
                    anchorEndIndex = i;
                    partnerStartIndex = j;
                    partnerEndIndex = j;
                    break;
                }
            }
            
            if (anchorStartIndex > 0) {
                break;
            }
        }
        
        if (anchorStartIndex > 0) {
            int midPointIndex = (anchorStartIndex + partnerStartIndex)/2;
            for (int i=anchorStartIndex-1; i>midPointIndex; --i) {
                Kmer seed = seqKmers.get(i);
                long rcHashVal = seed.getReverseComplementHash();

                for (int j=partnerEndIndex+1; j<=midPointIndex; ++j) {
                    Kmer candidate = seqKmers.get(j);
                    if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                        anchorStartIndex = i;
                        partnerEndIndex = j;
                        break;
                    }
                }
            }

            if (anchorStartIndex > 0 && anchorEndIndex-anchorStartIndex >= minMatchLen && partnerEndIndex-partnerStartIndex >= minMatchLen) {
                ++anchorEndIndex;
                ++partnerEndIndex;

                if (stranded) {
                    float partnerCov = getMedianKmerCoverage(seqKmers, partnerStartIndex, partnerEndIndex);
                    float middleCov = partnerEndIndex < anchorStartIndex ? getMedianKmerCoverage(seqKmers, partnerEndIndex, anchorStartIndex) : 0;
                    float anchorCov = getMedianKmerCoverage(seqKmers, anchorStartIndex, anchorEndIndex);

                    if (partnerCov > anchorCov) {
                        if (middleCov > anchorCov && middleCov >= partnerCov * maxCovGradient) {
                            return new ArrayList<>(seqKmers.subList(0, anchorStartIndex));
                        }
                        else {
                            return new ArrayList<>(seqKmers.subList(0, partnerEndIndex));
                        }
                    }
                    else {
                        if (middleCov > partnerCov && middleCov >= anchorCov * maxCovGradient) {
                            return new ArrayList<>(seqKmers.subList(partnerEndIndex, numKmers));
                        }
                        else {
                            return new ArrayList<>(seqKmers.subList(anchorStartIndex, numKmers));
                        }
                    }
                }
                else {
                    return new ArrayList<>(seqKmers.subList(partnerEndIndex,anchorStartIndex));
                }
            }
        }
        
        return seqKmers;
    }
    
    public static ArrayList<Kmer> trimReverseComplementArtifact(ArrayList<Kmer> seqKmers,
            int maxEdgeClip, int maxIndelSize, float minPercentIdentity, BloomFilterDeBruijnGraph graph) {
        int k = graph.getK();
        
        // search from left to right
        int numKmers = seqKmers.size();
        
        int leftIndex = -1;
        int rightIndex = -1;
        int indexToStop = Math.min(maxEdgeClip, numKmers);
        for (int i=0; i<indexToStop; ++i) {
            Kmer seed = seqKmers.get(i);
            long rcHashVal = seed.getReverseComplementHash();
            for (int j=i+1; j<numKmers; ++j) {
                Kmer candidate = seqKmers.get(j);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    leftIndex = i;
                    rightIndex = j;
                    break;
                }
            }
            
            if (leftIndex > 0) {
                break;
            }
        }
        
        if (leftIndex > 0 && rightIndex-leftIndex >= k) {
            int cutIndex = leftIndex+1;
            for (int i=k; i<rightIndex-leftIndex; i+=k) {
                Kmer seed = seqKmers.get(leftIndex+i);
                long rcHashVal = seed.getReverseComplementHash();
                Kmer candidate = seqKmers.get(rightIndex-i);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    cutIndex = leftIndex + i;
                }
                else {
                    break;
                }
            }
            
            for (int i=cutIndex-leftIndex; i<rightIndex-leftIndex; ++i) {
                Kmer seed = seqKmers.get(leftIndex+i);
                long rcHashVal = seed.getReverseComplementHash();
                Kmer candidate = seqKmers.get(rightIndex-i);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    cutIndex = leftIndex + i;
                }
                else {
                    break;
                }
            }
            
            cutIndex = Math.min(cutIndex, (leftIndex + rightIndex)/2);
            
            if (rightIndex >= numKmers-maxEdgeClip) {
                int cutLength = cutIndex - leftIndex;
                float leftMinCov = getMinimumKmerCoverage(seqKmers, cutIndex, Math.min(numKmers, cutIndex+k));
                float rightMinCov = getMinimumKmerCoverage(seqKmers, Math.max(0, numKmers-cutLength-k), numKmers-cutLength);
                if (leftMinCov > rightMinCov) {
                    seqKmers = new ArrayList<>(seqKmers.subList(0, Math.max(1, numKmers-cutLength-k)));
                }
                else {
                    seqKmers = new ArrayList<>(seqKmers.subList(cutIndex, numKmers));
                }
            }
            else {
                seqKmers = new ArrayList<>(seqKmers.subList(Math.min(numKmers-1, cutIndex+k), numKmers));
            }
        }
        
        // search from right to left
        numKmers = seqKmers.size();
        
        leftIndex = -1;
        rightIndex = -1;
        indexToStop = Math.max(0, numKmers-1-maxEdgeClip);
        for (int i=numKmers-1; i>=indexToStop; --i) {
            Kmer seed = seqKmers.get(i);
            long rcHashVal = seed.getReverseComplementHash();
            for (int j=i-1; j>=0; --j) {
                Kmer candidate = seqKmers.get(j);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    rightIndex = i;
                    leftIndex = j;
                    break;
                }
            }
            
            if (rightIndex > 0) {
                break;
            }
        }
        
        if (rightIndex > 0 && rightIndex-leftIndex >= k) {
            int cutIndex = rightIndex-1;
            for (int i=k; i<rightIndex-leftIndex; i+=k) {
                Kmer seed = seqKmers.get(rightIndex-i);
                long rcHashVal = seed.getReverseComplementHash();
                Kmer candidate = seqKmers.get(leftIndex+i);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    cutIndex = rightIndex - i;
                }
                else {
                    break;
                }
            }
            
            for (int i=rightIndex-cutIndex; i<rightIndex-leftIndex; ++i) {
                Kmer seed = seqKmers.get(rightIndex-i);
                long rcHashVal = seed.getReverseComplementHash();
                Kmer candidate = seqKmers.get(leftIndex+i);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    cutIndex = rightIndex - i;
                }
                else {
                    break;
                }
            }
            
            cutIndex = Math.max(cutIndex, (rightIndex + leftIndex)/2);
            
            if (leftIndex < maxEdgeClip) {
                int cutLength = rightIndex - cutIndex;
                float leftMinCov = getMinimumKmerCoverage(seqKmers, cutLength, Math.min(cutLength+k, numKmers));
                float rightMinCov = getMinimumKmerCoverage(seqKmers, Math.max(0, cutIndex-k), cutIndex);
                if (leftMinCov > rightMinCov) {
                    seqKmers = new ArrayList<>(seqKmers.subList(0, Math.max(1, cutIndex-k)));
                }
                else {
                    seqKmers = new ArrayList<>(seqKmers.subList(Math.min(cutLength+k, numKmers-1), numKmers));
                }
            }
            else {
                seqKmers = new ArrayList<>(seqKmers.subList(0, Math.max(1, cutIndex-k)));
            }
        }
        
        return seqKmers;
    }
    
    public static ArrayList<Kmer> trimHairpinBySequenceMatching(ArrayList<Kmer> seqKmers, int k, float minPercentIdentity, BloomFilterDeBruijnGraph graph) {
        int numKmers = seqKmers.size();
        int halfNumKmers = numKmers/2;
        
        final int maxSeedSearchDepth = Math.min(halfNumKmers, 200);
        final int maxLoopLength = Math.max(200, halfNumKmers);
        final int maxLoopDiameter = maxLoopLength/2;
        
        int lastSeedIndex = maxSeedSearchDepth;
        
        for (int i=0; i<lastSeedIndex; i+=k) {
            Kmer seed = seqKmers.get(i);
            long rcHashVal = seed.getReverseComplementHash();
            int rcIndex = -1;
            
            for (int j=i+1; j<numKmers; ++j) {
                Kmer candidate = seqKmers.get(j);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    rcIndex = j;
                    break;
                }
            }
            
            if (rcIndex >= 0) {
                int endIndex = rcIndex;
                int halfIndex = (i + endIndex)/2;
                
                if (i >= rcIndex-maxLoopLength) {
                    if (halfIndex < halfNumKmers) {
                        return new ArrayList<>(seqKmers.subList(halfIndex, numKmers));
                    }
                    else {
                        return new ArrayList<>(seqKmers.subList(0, halfIndex));
                    }
                }
                else {
                    int testLength = halfIndex-maxLoopDiameter+1-i;
                    byte[] left = graph.assembleBytes(seqKmers, i, halfIndex-maxLoopDiameter+1);
                    byte[] right = graph.assembleReverseComplementBytes(seqKmers, endIndex+1-testLength, endIndex+1);
                    
                    float pid = getPercentIdentity(left, right);
                    
                    if (pid >= minPercentIdentity) {
                        if (halfIndex < halfNumKmers) {
                            return new ArrayList<>(seqKmers.subList(halfIndex, numKmers));
                        }
                        else {
                            return new ArrayList<>(seqKmers.subList(0, halfIndex));
                        }
                    }
                }
                
                break;
            }
        }
        
        lastSeedIndex = numKmers-maxSeedSearchDepth;
        
        for (int i=numKmers-1; i>=lastSeedIndex; i-=k) {
            Kmer seed = seqKmers.get(i);
            long rcHashVal = seed.getReverseComplementHash();
            int rcIndex = -1;
            
            for (int j=i-1; j>=0; --j) {
                Kmer candidate = seqKmers.get(j);
                if (candidate.getHash() == rcHashVal && isReverseComplement(seed.bytes, candidate.bytes)) {
                    rcIndex = j;
                    break;
                }
            }

            if (rcIndex >= 0 && rcIndex < i) {
                int halfIndex = (rcIndex + i)/2;
                
                if (rcIndex >= i-maxLoopLength) {
                    if (halfIndex < halfNumKmers) {
                        return new ArrayList<>(seqKmers.subList(halfIndex, numKmers));
                    }
                    else {
                        return new ArrayList<>(seqKmers.subList(0, halfIndex));
                    }
                }
                else {
                    int testLength = halfIndex-maxLoopDiameter-rcIndex;

                    byte[] left = graph.assembleBytes(seqKmers, rcIndex, halfIndex-maxLoopDiameter+1);
                    byte[] right = graph.assembleReverseComplementBytes(seqKmers, i-testLength, i+1);
                    float pid = getPercentIdentity(left, right);
                    
                    if (pid >= minPercentIdentity) {
                        if (halfIndex < halfNumKmers) {
                            return new ArrayList<>(seqKmers.subList(halfIndex, numKmers));
                        }
                        else {
                            return new ArrayList<>(seqKmers.subList(0, halfIndex));
                        }
                    }
                }
                
                break;
            }
        }
        
        return seqKmers;
    }
    
    public static ArrayList<Kmer> trimHairpinByCoverage(ArrayList<Kmer> seqKmers, int k, float minPercentIdentity) {                
        ArrayDeque<Integer> minValIndexes = getMinCovIndexes(seqKmers, k);
                
        while (!minValIndexes.isEmpty()) {
            int numKmers = seqKmers.size();
            int halfIndex = numKmers/2;
            
            int clippingIndex = getReverseComplementArtifactClippingIndex(seqKmers, k, minPercentIdentity, minValIndexes);
            
            if (clippingIndex < 0) {
                break;
            }
            else if (clippingIndex > halfIndex) {
                seqKmers = new ArrayList<>(seqKmers.subList(0, clippingIndex));
            }
            else {
                seqKmers = new ArrayList<>(seqKmers.subList(clippingIndex, numKmers));
            }
            
            if (minValIndexes.size() == 1) {
                break;
            }
            
            minValIndexes = getMinCovIndexes(seqKmers, k);
        }
        
        return seqKmers;
    }
        
    private static ArrayDeque<Integer> getMinCovIndexes(ArrayList<Kmer> seqKmers, int k) {
        int numKmers = seqKmers.size();
        
        float min = Float.MAX_VALUE;
        ArrayDeque<Integer> minValIndexes = new ArrayDeque<>();
        
        for (int i=0; i<numKmers; ++i) {
            float c = seqKmers.get(i).count;
            if (c < min) {
                min = c;
                minValIndexes = new ArrayDeque<>();
                minValIndexes.add(i);
            }
            else if (c == min) {
                if (minValIndexes.getLast() < i-k) {
                    minValIndexes.add(i);
                }
            }
        }
        
        if (!minValIndexes.isEmpty()) {
            if (minValIndexes.getFirst() == 0) {
                minValIndexes.removeFirst();
            }
        }
        
        if (!minValIndexes.isEmpty()) {
            if (minValIndexes.getLast() == numKmers-1) {
                minValIndexes.removeLast();
            }
        }
        
        return minValIndexes;
    }
    
    public static int getReverseComplementArtifactClippingIndex(ArrayList<Kmer> seqKmers, int k, float minPercentIdentity, ArrayDeque<Integer> minValIndexes) {        
        int numKmers = seqKmers.size();
        int halfIndex = numKmers/2;

        if (minValIndexes != null && !minValIndexes.isEmpty()) {
            float min = seqKmers.get(minValIndexes.getFirst()).count;
            
            for (int bendIndex : minValIndexes) {
                if (bendIndex < halfIndex) {
                    int numNeeded = (int) Math.floor(minPercentIdentity*(bendIndex+1));
                    if (numNeeded >= 1) {
                        for (int i=bendIndex; i<numKmers; ++i) {
                            Kmer kmer = seqKmers.get(i);
                            if (kmer.count == min) {
                                bendIndex = i;
                            }
                            else {
                                break;
                            }
                        }
                        
                        HashSet<Long> hashVals = new HashSet<>(numKmers - bendIndex + 1);
                        for (int i=bendIndex+1; i<numKmers; ++i) {
                            hashVals.add(seqKmers.get(i).getHash());
                        }
                        
                        for (int i=0; i<=bendIndex; ++i) {
                            if (hashVals.contains(seqKmers.get(i).getReverseComplementHash())) {
                                if (--numNeeded <= 0) {
                                    break;
                                }
                            }
                        }

                        if (numNeeded <= 0) {
                            return bendIndex;
                        }
                    }
                }
                else {
                    int numNeeded = (int) Math.floor(minPercentIdentity*(numKmers-bendIndex+1));
                    if (numNeeded >= 1) {
                        for (int i=bendIndex-1; i>=0; --i) {
                            Kmer kmer = seqKmers.get(i);
                            if (kmer.count == min) {
                                bendIndex = i;
                            }
                            else {
                                break;
                            }
                        }
                        
                        HashSet<Long> hashVals = new HashSet<>(bendIndex + 1);
                        for (int i=0; i<bendIndex; ++i) {
                            hashVals.add(seqKmers.get(i).getHash());
                        }
                        
                        for (int i=bendIndex; i<numKmers; ++i) {
                            if (hashVals.contains(seqKmers.get(i).getReverseComplementHash())) {
                                if (--numNeeded <= 0) {
                                    break;
                                }
                            }
                        }

                        if (numNeeded <= 0) {
                            return bendIndex;
                        }
                    }
                }
            }
        }
        
        return -1;
    }
        
    public static boolean isTemplateSwitch2(ArrayList<Kmer> seqKmers, BloomFilterDeBruijnGraph graph, BloomFilter assembledKmers, int lookahead, float minPercentIdentity) {
        int numKmers = seqKmers.size();
        int k = graph.getK();
        int maxNumKmersInLoop = 2*k;
        
        if (assembledKmers.lookup(seqKmers.get(numKmers-1).getHash())) {
            int start = Math.max(0, numKmers-2);
            for (; start>=0; --start) {
                if (!assembledKmers.lookup(seqKmers.get(start).getHash())) {
                    if (start-k >=0) {
                        ArrayDeque<Kmer> path = getMaxCoveragePath(graph, 
                                                    seqKmers.get(start-k), 
                                                    seqKmers.get(start+1),
                                                    k+1, 
                                                    lookahead, 
                                                    assembledKmers);
                        if (path != null) {
                            start -= k; 
                            continue;
                        }
                    }
                    
                    ++start;
                    break;
                }
            }
            
            if (start < k) {
                // no unassembled kmers
                return false;
            }
            
            float medCoverage = getMedianKmerCoverage(seqKmers, start, numKmers);

            int backBoneKmerIndex = -1;
            for (int i=start; i<numKmers; ++i) {
                if (seqKmers.get(i).count >= medCoverage) {
                    backBoneKmerIndex = i;
                    break;
                }
            }

            if (backBoneKmerIndex >= 0) {
                Kmer backBoneKmer = seqKmers.get(backBoneKmerIndex);

                ArrayDeque<Kmer> extension = greedyExtendLeft(graph, backBoneKmer, lookahead, 1000, assembledKmers);
                extension.add(backBoneKmer);
                extension.addAll(greedyExtendRight(graph, backBoneKmer, lookahead, 1000, assembledKmers));
                
                String backBone = graph.assemble(extension);

                String tipRC = reverseComplement(graph.assemble(seqKmers, 0, Math.max(1, start-maxNumKmersInLoop)));

                if (backBone.contains(tipRC)) {
                    return true;
                }

//                if (found) {
//                    System.out.println(">seq\n" + graph.assemble(seqKmers));
//                    System.out.println(">backbone\n" + backBone);
//                    System.out.println(">tipRC\n" + tipRC);
//                    System.out.println("yoohoo");
//                }
            }
        }
        
        if (assembledKmers.lookup(seqKmers.get(0).getHash())) {
            int end = Math.min(1, numKmers-1);
            for (; end<numKmers; ++end) {
                if (!assembledKmers.lookup(seqKmers.get(end).getHash())) {
                    if (end+k < numKmers) {
                        ArrayDeque<Kmer> path = getMaxCoveragePath(graph, 
                                                    seqKmers.get(end-1), 
                                                    seqKmers.get(end+k),
                                                    k+1, 
                                                    lookahead, 
                                                    assembledKmers);
                        if (path != null) {
                            end += k;
                            continue;
                        }
                    }
                    
                    break;
                }
            }
            
            if (end >= numKmers-k) {
                // no unassembled kmers
                return false;
            }
            
            float medCoverage = getMedianKmerCoverage(seqKmers, 0, end);

            int backBoneKmerIndex = -1;
            for (int i=0; i<end; ++i) {
                if (seqKmers.get(i).count >= medCoverage) {
                    backBoneKmerIndex = i;
                    break;
                }
            }

            if (backBoneKmerIndex >= 0) {
                Kmer backBoneKmer = seqKmers.get(backBoneKmerIndex);

                ArrayDeque<Kmer> extension = greedyExtendLeft(graph, backBoneKmer, lookahead, 1000, assembledKmers);
                extension.add(backBoneKmer);
                extension.addAll(greedyExtendRight(graph, backBoneKmer, lookahead, 1000, assembledKmers));                

                String backBone = graph.assemble(extension);

                String tipRC = reverseComplement(graph.assemble(seqKmers, Math.min(numKmers-1, end+maxNumKmersInLoop), numKmers));

                if (backBone.contains(tipRC)) {
                    return true;
                }

//                if (found) {
//                    System.out.println(">seq\n" + graph.assemble(seqKmers));
//                    System.out.println(">backbone\n" + backBone);
//                    System.out.println(">tipRC\n" + tipRC);
//                    System.out.println("yoohoo");
//                }
            }
        }
        
        return false;
    }
    
    public static boolean isTemplateSwitch(ArrayList<Kmer> seqKmers, BloomFilterDeBruijnGraph graph, BloomFilter assembledKmers, int lookahead) {
        int k = graph.getK();   
        int numKmers = seqKmers.size();
        
        float leftEdgeCov = getMinimumKmerCoverage(seqKmers, 0, Math.min(lookahead, numKmers));
        float rightEdgeCov = getMinimumKmerCoverage(seqKmers, Math.max(0, numKmers-lookahead), numKmers);
        
        if (assembledKmers.lookup(seqKmers.get(0).getHash()) &&
                (!assembledKmers.lookup(seqKmers.get(numKmers-1).getHash()) || leftEdgeCov > rightEdgeCov)) {
            int i = 1;
            for (; i < numKmers; ++i) {
                if (!assembledKmers.lookup(seqKmers.get(i).getHash())) {
                    break;
                }
            }
            
            if (i == numKmers) {
                return false;
            }
            
            --i;
            
            String tipRC = reverseComplement(graph.assemble(seqKmers, Math.min(i+k, numKmers-1), numKmers));

            if (i-lookahead >= 0) {
                i -= lookahead;
            }

            ArrayDeque<Kmer> leftExtension = greedyExtendLeft(graph, seqKmers.get(0), lookahead, 1000, assembledKmers);
            ArrayDeque<Kmer> rightExtension = greedyExtendRight(graph, seqKmers.get(i), lookahead, 1000, assembledKmers);

            leftExtension.addAll(seqKmers.subList(0, i+1));
            leftExtension.addAll(rightExtension);

            String backbone = graph.assemble(leftExtension);
            if (backbone.contains(tipRC)) {
                return true;
            }
//            else {
//                HashSet<Kmer> intersection = new HashSet<>(rightExtension);
//                intersection.retainAll(greedyExtendLeft(graph, graph.getKmer(getFirstKmer(tipRC, k)), lookahead, 1000, assembledKmers));
//                if (!intersection.isEmpty()) {
//                    return true;
//                }
//
//                intersection = new HashSet<>(leftExtension);
//                intersection.retainAll(greedyExtendRight(graph, graph.getKmer(getLastKmer(tipRC, k)), lookahead, 1000, assembledKmers));
//                if (!intersection.isEmpty()) {
//                    return true;
//                }                
//            }
        }
        else if (assembledKmers.lookup(seqKmers.get(numKmers-1).getHash()) &&
                (!assembledKmers.lookup(seqKmers.get(0).getHash()) || leftEdgeCov < rightEdgeCov)) {
            int j = numKmers-2;
            for (; j >= 0; --j) {
                if (!assembledKmers.lookup(seqKmers.get(j).getHash())) {
                    break;
                }
            }
            
            if (j == -1) {
                return false;
            }
            
            ++j;
            
            String tipRC = reverseComplement(graph.assemble(seqKmers, 0, Math.max(1, j-k)));

            if (j+lookahead < numKmers) {
                j += lookahead;
            }

            ArrayDeque<Kmer> leftExtension = greedyExtendLeft(graph, seqKmers.get(j), lookahead, 1000, assembledKmers);                
            ArrayDeque<Kmer> rightExtension = greedyExtendRight(graph, seqKmers.get(numKmers-1), lookahead, 1000, assembledKmers);

            leftExtension.addAll(seqKmers.subList(j, numKmers));
            leftExtension.addAll(rightExtension);

            String backbone = graph.assemble(leftExtension);
            if (backbone.contains(tipRC)) {
                return true;
            }
//            else {
//                HashSet<Kmer> intersection = new HashSet<>(rightExtension);
//                intersection.retainAll(greedyExtendLeft(graph, graph.getKmer(getFirstKmer(tipRC, k)), lookahead, 1000, assembledKmers));
//                if (!intersection.isEmpty()) {
//                    return true;
//                }
//
//                intersection = new HashSet<>(leftExtension);
//                intersection.retainAll(greedyExtendRight(graph, graph.getKmer(getLastKmer(tipRC, k)), lookahead, 1000, assembledKmers));
//                if (!intersection.isEmpty()) {
//                    return true;
//                }                
//            }
        }
        
        return false;
    }
    
    public static boolean isBluntEndArtifact(ArrayList<Kmer> seqKmers, BloomFilterDeBruijnGraph graph, BloomFilter assembledKmers, int maxDepth) {
        if (maxDepth <= 0) {
            return false;
        }
        
        int numKmers = seqKmers.size();
        int d = graph.getReadPairedKmerDistance();
        
        float leftEdgeCov = getMinimumKmerCoverage(seqKmers, 0, Math.min(maxDepth, numKmers));
        float rightEdgeCov = getMinimumKmerCoverage(seqKmers, Math.max(0, numKmers-maxDepth), numKmers);
        
        if (assembledKmers.lookup(seqKmers.get(0).getHash()) &&
                (!assembledKmers.lookup(seqKmers.get(numKmers-1).getHash()) || leftEdgeCov > rightEdgeCov)) {
            int i = 1;
            for (; i < numKmers; ++i) {
                if (!assembledKmers.lookup(seqKmers.get(i).getHash())) {
                    break;
                }
            }
            
            if (i == numKmers || i < numKmers-d) {
                return false;
            }

            if (!hasDepthRight(seqKmers.get(numKmers-1), graph, maxDepth) &&
                    getMedianKmerCoverage(seqKmers, 0, i) > getMedianKmerCoverage(seqKmers, i, numKmers) &&
                    hasDepthRight(seqKmers.get(i-1), graph, numKmers-i, assembledKmers)) {
                return true;
            }
        }
        else if (assembledKmers.lookup(seqKmers.get(numKmers-1).getHash()) &&
                (!assembledKmers.lookup(seqKmers.get(0).getHash()) || leftEdgeCov < rightEdgeCov)) {
            int j = numKmers-2;
            for (; j >= 0; --j) {
                if (!assembledKmers.lookup(seqKmers.get(j).getHash())) {
                    break;
                }
            }
            
            if (j == -1 || j > d) {
                return false;
            }
            
            if (!hasDepthLeft(seqKmers.get(0), graph, maxDepth) &&
                    getMedianKmerCoverage(seqKmers, j+1, numKmers) > getMedianKmerCoverage(seqKmers, 0, j+1) &&
                    hasDepthLeft(seqKmers.get(j+1), graph, j+1, assembledKmers)) {
                return true;
            }
        }
        
        return false;
    }
    
//    public static void main(String[] args) {
////        String seq = "AAAAAAAAAAA";
//        String seq = "AAAAAAAAAAACCC";
////        String seq = "AAAAAAAAAAACCCCCCCCCCCGGGGGGGGGGG";
////        String seq = "AAAAAAAAAAACCCCCCCCCCCGGGGGGGGGGGTTT";
//        System.out.println(seq);
//        
//        int k = 11;
//        
//        ArrayList<Kmer2> kmers = new ArrayList<>();
//        for (String kmer : SeqUtils.kmerize(seq, k)) {
//             kmers.add(new Kmer(kmer, 1, new long[0]));
//        }
//
//        System.out.println(assemble(kmers, k));
//    }
}
