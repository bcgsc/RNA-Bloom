/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import rnabloom.RNABloom;
import rnabloom.RNABloom.ReadPair;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
import static rnabloom.util.SeqUtils.getFirstKmer;
import static rnabloom.util.SeqUtils.getLastKmer;
import static rnabloom.util.SeqUtils.getPercentIdentity;
import static rnabloom.util.SeqUtils.isLowComplexity2;
import static rnabloom.util.SeqUtils.overlapMaximally;
import static rnabloom.util.SeqUtils.stringToBytes;
import static rnabloom.util.SeqUtils.isLowComplexity2;
import static rnabloom.util.SeqUtils.isLowComplexity2;
import static rnabloom.util.SeqUtils.isLowComplexity2;

/**
 *
 * @author gengar
 */
public final class GraphUtils {

    private static float getMinimumKmerCoverage(final Iterable<Kmer> kmers) {
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
    
    private static float getMinimumKmerCoverage(final ArrayList<Kmer> kmers, int start, int end) {
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
    
    public static float getMedianKmerCoverage(final ArrayList<Kmer> kmers, int start, int end) {
        int range = end-start;
        
        float[] covs = new float[range];
        for (int i=0; i<range; ++i) {
            covs[i] = kmers.get(start+i).count;
        }

        return getMedian(covs);
    }
    
    public static float getMedianKmerCoverage(final ArrayDeque<Kmer> kmers) {
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
        
//        ArrayList<Float> counts = new ArrayList<>(numKmers);
//        for (Kmer kmer : kmers) {
//            counts.add(kmer.count);
//        }
//        
//        Collections.sort(counts);
//        
//        if (numKmers % 2 == 0) {
//            return (counts.get(halfNumKmers) + counts.get(halfNumKmers -1))/2.0f;
//        }
//        
//        return counts.get(halfNumKmers);
    }
    
    public static float getMaxMedianCoverageRight(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead) {
        ArrayDeque<Kmer> neighbors = graph.getSuccessors(source);
        
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
                    neighbors = graph.getSuccessors(cursor);
                    if (!neighbors.isEmpty()) {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        frontier.add(neighbors);
                        continue;
                    }
                }

                if (path.size() == lookahead) {
                    // we only calculate coverage if path is long enough
                    float pathCov = getMedianKmerCoverage(path);
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
        ArrayDeque<Kmer> neighbors = graph.getSuccessors(source, bf);
        
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
                    neighbors = graph.getSuccessors(cursor, bf);
                    if (!neighbors.isEmpty()) {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        frontier.add(neighbors);
                        continue;
                    }
                }

                if (path.size() == lookahead) {
                    // we only calculate coverage if path is long enough
                    float pathCov = getMedianKmerCoverage(path);
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
        ArrayDeque<Kmer> neighbors = graph.getPredecessors(source);
        
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
                    neighbors = graph.getPredecessors(cursor);
                    if (!neighbors.isEmpty()) {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        frontier.add(neighbors);
                        continue;
                    }
                }

                if (path.size() == lookahead) {
                    // we only calculate coverage if path is long enough
                    float pathCov = getMedianKmerCoverage(path);
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
        ArrayDeque<Kmer> neighbors = graph.getPredecessors(source, bf);
        
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
                    neighbors = graph.getPredecessors(cursor, bf);
                    if (!neighbors.isEmpty()) {
                        cursor = neighbors.removeFirst();
                        path.add(cursor);
                        frontier.add(neighbors);
                        continue;
                    }
                }

                if (path.size() == lookahead) {
                    // we only calculate coverage if path is long enough
                    float pathCov = getMedianKmerCoverage(path);
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
                }
                return bestKmer;
            }
        }
    }
    
    public static Kmer greedyExtendRightOnce(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead) {
        return greedyExtendRightOnce(graph, graph.getSuccessors(source), lookahead);
    }
    
    public static Kmer greedyExtendRightOnce(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead, BloomFilter bf) {
        return greedyExtendRightOnce(graph, graph.getSuccessors(source), lookahead, bf);
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
                }
                return bestKmer;
            }
        }
    }
    
    public static Kmer greedyExtendLeftOnce(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead) {
        return greedyExtendLeftOnce(graph, graph.getPredecessors(source), lookahead);
    }

    public static Kmer greedyExtendLeftOnce(final BloomFilterDeBruijnGraph graph, final Kmer source, final int lookahead, BloomFilter bf) {
        return greedyExtendLeftOnce(graph, graph.getPredecessors(source), lookahead, bf);
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
                }
                return bestKmer;
            }
        }
    }
    
    public static boolean represented(final ArrayList<Kmer> kmers,
                                    final BloomFilterDeBruijnGraph graph,
                                    final BloomFilter bf,
                                    final int lookahead,
                                    final int maxIndelSize,
                                    final float percentIdentity) {
        int numKmers = kmers.size();
        int maxIndex = numKmers - 1;
        
        int k = graph.getK();
        int maxNumBubbleKmers = 3*k;
        
        int lastRepresentedKmerFoundIndex = -1;
        
        for (int i=0; i<numKmers; ++i) {
            
            if (bf.lookup(kmers.get(i).hashVals)) {
                int startIndex = i;
                int endIndex = i;
                for (int j=i+1; j<numKmers; ++j) {
                    if (bf.lookup(kmers.get(j).hashVals)) {
                        endIndex = j;
                    }
                    else {
                        break;
                    }
                }
                
                int assembledRange = endIndex - startIndex + 1;
                
                if (assembledRange >= 3) {                    
                    if (startIndex > 0) {
                                                
                        if (lastRepresentedKmerFoundIndex < 0) {                            
                            if (startIndex >= k) {
                                // check left edge kmers
                                ArrayDeque testEdgeKmers = greedyExtendLeft(graph, kmers.get(startIndex), lookahead, startIndex, bf);
                                if (testEdgeKmers.size() != startIndex ||
                                        getPercentIdentity(assemble(testEdgeKmers, k), assemble(kmers, k, 0, startIndex)) < percentIdentity) {
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
                                    if (left == 0 || !bf.lookup(kmers.get(--left).hashVals)) {
                                        break;
                                    }
                                }

                                for (int j=0; j<numMissing; ++j) {
                                    if (right == maxIndex || !bf.lookup(kmers.get(++right).hashVals)) {
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
                                        getPercentIdentity(assemble(testPathKmers, k), assemble(kmers, k, left+1, right)) < percentIdentity)) {
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
                if (expectedLen >= k) {                    
                    ArrayDeque testEdgeKmers = greedyExtendRight(graph, kmers.get(lastRepresentedKmerFoundIndex), lookahead, expectedLen, bf);
                    if (testEdgeKmers.size() != expectedLen ||
                            getPercentIdentity(assemble(testEdgeKmers, k), assemble(kmers, k, lastRepresentedKmerFoundIndex+1, numKmers)) < percentIdentity) {
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
        
        ArrayDeque<Kmer> frontier = new ArrayDeque<>();
        frontier.addAll(graph.getSuccessors(left));
        
        HashSet<String> kmersInFrontier = new HashSet<>();
        ArrayDeque<Kmer> newFrontier;
        for (int i=1; i<lowerBound; ++i) {
            kmersInFrontier.clear();
            newFrontier = new ArrayDeque<>();
            for (Kmer kmer : frontier) {
                if (bf.lookup(kmer.hashVals)) {
                    for (Kmer s : graph.getSuccessors(kmer)) {
                        String seq = s.toString();
                        if (!kmersInFrontier.contains(seq)) { 
                            newFrontier.add(s);
                            kmersInFrontier.add(seq);
                        }
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
                if (bf.lookup(kmer.hashVals)) {
                    if (kmer.equals(right)) {
                        return true;
                    }
                    newFrontier.add(kmer);
                    for (Kmer s : graph.getSuccessors(kmer)) {
                        String seq = s.toString();
                        if (!kmersInFrontier.contains(seq)) { 
                            newFrontier.add(s);
                            kmersInFrontier.add(seq);
                        }
                    }
                }
            }
            
            if (newFrontier.isEmpty()) {
                return false;
            }
            
            frontier = newFrontier;
        }
        
        for (Kmer kmer : frontier) {
            if (bf.lookup(kmer.hashVals)) {
                if (kmer.equals(right)) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    /**
     * 
     * @param graph
     * @param left
     * @param right
     * @param bound
     * @param lookahead
     * @return 
     */
    public static ArrayDeque<Kmer> getMaxCoveragePath(BloomFilterDeBruijnGraph graph, Kmer left, Kmer right, int bound, int lookahead) {
        
        HashSet<String> leftPathKmers = new HashSet<>(bound);
        
        /* extend right */
        ArrayDeque<Kmer> leftPath = new ArrayDeque<>(bound);
        Kmer best;
        ArrayDeque<Kmer> neighbors;
        
        best = left;

        for (int depth=0; depth < bound; ++depth) {
            neighbors = graph.getSuccessors(best);
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
                    String seq = best.toString();
                    if (leftPathKmers.contains(seq)) {
                        break;
                    }
                    else {
                        leftPathKmers.add(seq);
                        leftPath.add(best);
                    }
                }
            }
        }
        
        HashSet<String> rightPathKmers = new HashSet<>(bound);
        
        /* not connected, search from right */
        ArrayDeque<Kmer> rightPath = new ArrayDeque<>(bound);
        best = right;
        for (int depth=0; depth < bound; ++depth) {
            neighbors = graph.getPredecessors(best);
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
                    String bestSeq = best.toString();
                    if (rightPathKmers.contains(bestSeq)) {
                        return null;
                    }
                    else if (leftPathKmers.contains(bestSeq)) {
                        /* right path intersects the left path */
                        
                        if (isLowComplexity2(bestSeq)) {
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
                        rightPathKmers.add(bestSeq);
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
               
        HashSet<String> leftPathKmers = new HashSet<>(bound);
        
        /* extend right */
        ArrayDeque<Kmer> leftPath = new ArrayDeque<>(bound);
        Kmer best;
        ArrayDeque<Kmer> neighbors;
        
        best = left;

        for (int depth=0; depth < bound; ++depth) {
            neighbors = graph.getSuccessors(best, bf);
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
                    String bestSeq = best.toString();
                    
                    if (leftPathKmers.contains(bestSeq)) {
                        break;
                    }
                    else {
                        leftPathKmers.add(bestSeq);
                        leftPath.add(best);
                    }
                }
            }
        }
        
        HashSet<String> rightPathKmers = new HashSet<>(bound);
        
        /* not connected, search from right */
        ArrayDeque<Kmer> rightPath = new ArrayDeque<>(bound);
        best = right;
        for (int depth=0; depth < bound; ++depth) {
            neighbors = graph.getPredecessors(best, bf);
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
                
                String bestSeq = best.toString();
                
                if (best.equals(left)) {
                    return rightPath;
                }
                else if (leftPathKmers.contains(bestSeq)) {
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
                else if (!rightPathKmers.contains(bestSeq)) {
                    rightPathKmers.add(bestSeq);
                    rightPath.addFirst(best);
                }
                else {
                    return null;
                }
            }
        }
        
        return null;
    }
    
    public static float getMedianModifyInput(float[] arr) {
        int len = arr.length;
        Arrays.sort(arr);
        int halfLen = len/2;
        if (len % 2 == 0) {
            return (arr[halfLen-1] + arr[halfLen])/2.0f;
        }
        
        return arr[halfLen];
    }    
    
    public static float getMedian(float[] arr) {
        int len = arr.length;
        float[] a = Arrays.copyOf(arr, len);
        Arrays.sort(a);
        int halfLen = len/2;
        if (len % 2 == 0) {
            return (a[halfLen-1] + a[halfLen])/2.0f;
        }
        
        return a[halfLen];
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
        
    private static class Stats {
        public float min;
        public float q1;
        public float median;
        public float q3;
        public float max;
    }
    
    private static Stats getStats(final List<Float> arr) {
        ArrayList<Float> a = new ArrayList<>(arr);
        Collections.sort(a);
        
        Stats stats = new Stats();
        
        int arrLen = a.size();
        int halfLen = arrLen/2;
        int q1Index = arrLen/4;
        int q3Index = halfLen+q1Index;
        
        stats.min = a.get(0);
        stats.max = a.get(arrLen-1);
        
        if (arrLen % 2 == 0) {
            stats.median = (a.get(halfLen-1) + a.get(halfLen))/2.0f;
        }
        else {
            stats.median = a.get(halfLen);
        }
        
        if (arrLen % 4 == 0) {
            stats.q1 = (a.get(q1Index-1) + a.get(q1Index))/2.0f;
            stats.q3 = (a.get(q3Index-1) + a.get(q3Index))/2.0f;
        }
        else {
            stats.q1 = a.get(q1Index);
            stats.q3 = a.get(q3Index);
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
        
        return getMedianModifyInput(covs);
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
        
        return getMedianModifyInput(covs);
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
    
    public static String correctMismatches(String seq, BloomFilterDeBruijnGraph graph, int lookahead, int mismatchesAllowed) {
        int numCorrected = 0;
        int seqLen = seq.length();
        int k = graph.getK();
        
        if (seqLen < k) {
            // no correction
            return seq;
        }
        
        int numKmers = seqLen-k+1;
        StringBuilder sb = new StringBuilder(seq);
        
        float bestCov, cov;
        String kmer, guide;
        ArrayDeque<String> variants;
        char bestBase = 'N';
        int guideEnd = -1;
        int guideLen = -1;
        
        // correct from start
        for (int i=0; i<numKmers; ++i) {
            int j = i+k;
            kmer = sb.substring(i, j);
            variants = graph.getRightVariants(kmer);
            if (!variants.isEmpty()) {
                guideEnd = Math.min(j+lookahead, seqLen);
                guide = sb.substring(j, guideEnd);
                guideLen = guideEnd - j;
                bestCov = 0;
                
                if (graph.contains(kmer)) {
                    bestCov = rightGuidedMedianCoverage(graph, kmer, guide, guideLen);
                }
                
                boolean corrected = false;
                for (String v : variants) {
                    cov = rightGuidedMedianCoverage(graph, v, guide, guideLen);
                    if (cov > bestCov) {
                        bestCov = cov;
                        bestBase = v.charAt(k-1);
                        corrected = true;
                    }
                }
                
                if (corrected) {
                    if (++numCorrected > mismatchesAllowed) {
                        // too many mismatches
                        return seq;
                    }
                    else {
                        sb.setCharAt(j-1, bestBase);
                    }
                    
                    //i += lookahead;
                }
            }
        }
        
        // correct from end
        for (int i=seqLen-k; i>=0; --i) {
            kmer = sb.substring(i, i+k);
            variants = graph.getLeftVariants(kmer);
            if (!variants.isEmpty()) {
                guideEnd = Math.max(0, i-lookahead);
                guideLen = i - guideEnd;
                guide = sb.substring(guideEnd, i);
                bestCov = 0;
                
                if (graph.contains(kmer)) {
                    bestCov = leftGuidedMedianCoverage(graph, kmer, guide, guideLen);
                }
                
                boolean corrected = false;
                for (String v : variants) {
                    cov = leftGuidedMedianCoverage(graph, v, guide, guideLen);
                    if (cov > bestCov) {
                        bestCov = cov;
                        bestBase = v.charAt(0);
                        corrected = true;
                    }
                }
                
                if (corrected) {                    
                    if (++numCorrected > mismatchesAllowed) {
                        // too many mismatches
                        return seq;
                    }
                    else {
                        sb.setCharAt(i, bestBase);
                    }
                    
                    //i -= lookahead;
                }
            }
        }
        
        if (numCorrected == 0) {
            return seq;
        }
        
        String seq2 = sb.toString();
        if (graph.isValidSeq(seq2)) {
            return seq2;
        }
        
        return seq;
    }
    
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
    
    private static ArrayList<Kmer> correctErrorHelper(ArrayList<Kmer> kmers,
                                                    BloomFilterDeBruijnGraph graph, 
                                                    int lookahead,
                                                    int maxIndelSize,
                                                    float covThreshold,
                                                    float percentIdentity) {
        
        
        boolean corrected = false;
        int numKmers = kmers.size();
        int lastKmerIndex = numKmers-1;
        
        int k = graph.getK();
        
//        int kMinus1 = graph.getKMinus1();

        ArrayList<Kmer> kmers2 = new ArrayList<>(numKmers + maxIndelSize);
        int numBadKmersSince = 0;
        Kmer kmer;
        for (int i=0; i<numKmers; ++i) {
            kmer = kmers.get(i);
            
            if (kmer.count >= covThreshold &&
                    (i==0 || kmers.get(i-1).count >= covThreshold) &&
                    (i==lastKmerIndex || kmers.get(i+1).count >= covThreshold)) {
                
                if (numBadKmersSince > 0) {
                    if (kmers2.isEmpty()) {
                        // check left edge
                        ArrayDeque<Kmer> leftVars = graph.getLeftVariants(kmers.get(i-1));
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

//                            // extract tip sequence
//                            StringBuilder sb = new StringBuilder(i);
//                            for (int j=0; j<i-1; ++j) {
//                                sb.append(kmers.get(j).seq.charAt(0));
//                            }
//                            String tip = sb.toString();
//
//                            // identify the best variant tip
//                            String best = null;
//                            float bestCov = tipMedCov;
//                            for (Kmer v : leftVars) {
//                                String s = tip + v.seq;
//                                float[] vc = graph.getMinMedianMaxKmerCoverage(s);
//                                if (vc[0] > 0 && vc[1] > bestCov) {
//                                    best = s;
//                                    bestCov = vc[1];
//                                }
//                            }
//
//                            if (best != null) {
//                                corrected = true;
//                                // replace with best variant tip
//                                kmers2.addAll(graph.getKmers(best));
//                            }
//                            else {
                                ArrayDeque<Kmer> greedyTipKmers = greedyExtendLeft(graph, kmer, lookahead, numBadKmersSince);
                                if (greedyTipKmers.size() == numBadKmersSince && getMedianKmerCoverage(greedyTipKmers) > tipMedCov) {
                                    if (getPercentIdentity(assemble(greedyTipKmers, k), assemble(kmers, k, 0, i)) >= percentIdentity){
                                        corrected = true;
                                        kmers2.addAll(greedyTipKmers);
                                    }
                                    else if (!graph.hasPredecessors(kmers.get(0)) && numBadKmersSince < k) {
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
//                            }
                        }
                    }
                    else {
                        ArrayDeque<Kmer> path = getMaxCoveragePath(graph, kmers2.get(kmers2.size()-1), kmer, numBadKmersSince + maxIndelSize, lookahead);
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
                                        getPercentIdentity(assemble(path, k), assemble(kmers, k, i-numBadKmersSince, i)) >= percentIdentity)) {
                                
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
            ArrayDeque<Kmer> rightVars = graph.getRightVariants(kmers.get(i));
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

//                // extract tip sequence
//                StringBuilder sb = new StringBuilder(numBadKmersSince);
//                for (int j=i+1; j<numKmers; ++j) {
//                    sb.append(kmers.get(j).seq.charAt(kMinus1));
//                }
//                String tip = sb.toString();
//
//                // identify the best variant tip
//                String best = null;
//                float bestCov = tipMedCov;
//                for (Kmer v : rightVars) {
//                    String s = v.seq + tip;
//                    float[] vc = graph.getMinMedianMaxKmerCoverage(s);
//                    if (vc[0] > 0 && vc[1] > bestCov) {
//                        best = s;
//                        bestCov = vc[1];
//                    }
//                }
//
//                if (best != null) {
//                    corrected = true;
//                    // replace with best variant tip
//                    kmers2.addAll(graph.getKmers(best));
//                }
//                else {
                    ArrayDeque<Kmer> greedyTipKmers = greedyExtendRight(graph, kmers.get(i-1), lookahead, numBadKmersSince);
                    if (greedyTipKmers.size() == numBadKmersSince && getMedianKmerCoverage(greedyTipKmers) > tipMedCov) {
                        if (getPercentIdentity(assemble(greedyTipKmers, k), assemble(kmers, k, i, numKmers)) >= percentIdentity){
                            corrected = true;
                            kmers2.addAll(greedyTipKmers);
                        }
                        else if (!graph.hasSuccessors(kmers.get(numKmers-1)) && numBadKmersSince < k) {
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
//                }
            }
        }
        
        if (corrected) {
            return kmers2;
        }
        
        // no changes
        return null;
    }
    
    public static ArrayList<Kmer> correctErrorsSE(ArrayList<Kmer> kmers,
                                    BloomFilterDeBruijnGraph graph, 
                                    int lookahead,
                                    int maxIndelSize,
                                    float maxCovGradient, 
                                    float covFPR,
                                    float percentIdentity) {
        
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
        
        float covThreshold = covs[startIndex];
        float c = -1;
        for (int i=startIndex-1; i>=0; --i) {
            c = covs[i];
            if (covThreshold * maxCovGradient > c) {
                thresholdFound = true;
                break;
            }
            covThreshold = c;
        }

        if (thresholdFound) {
            return correctErrorHelper(kmers,
                                        graph, 
                                        lookahead,
                                        maxIndelSize,
                                        covThreshold,
                                        percentIdentity);
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
                                            float percentIdentity) {
        
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
            int startIndex = numLeftKmers - 1 - numFalsePositivesAllowed;
            float leftCovThreshold = covs[startIndex];
            float c;
            for (int i=startIndex-1; i>=0; --i) {
                c = covs[i];
                if (leftCovThreshold * maxCovGradient > c) {
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
            startIndex = numRightKmers - 1 - numFalsePositivesAllowed;
            float rightCovThreshold = covs[startIndex];
            for (int i=startIndex-1; i>=0; --i) {
                c = covs[i];
                if (rightCovThreshold * maxCovGradient > c) {
                    rightThresholdFound = true;
                    break;
                }
                rightCovThreshold = c;
            }

            // kmers with coverage lower than theshold will be replaced
            if (leftThresholdFound || rightThresholdFound) {
                float covThreshold = Math.min(leftCovThreshold, rightCovThreshold);

                if (covThreshold >= minCovThreshold) {
                    // correct left read
                    ArrayList<Kmer> leftKmers2 = correctErrorHelper(leftKmers,
                                                                    graph, 
                                                                    lookahead,
                                                                    maxIndelSize,
                                                                    covThreshold,
                                                                    percentIdentity);

                    if (leftKmers2 != null) {
                        leftKmers = leftKmers2;
                        leftCorrected = true;
                    }

                    // correct right read
                    rightCorrected = false;
                    ArrayList<Kmer> rightKmers2 = correctErrorHelper(rightKmers,
                                                                    graph, 
                                                                    lookahead,
                                                                    maxIndelSize,
                                                                    covThreshold,
                                                                    percentIdentity);

                    if (rightKmers2 != null) {
                        rightKmers = rightKmers2;
                        rightCorrected = true;
                    }

                    if (leftKmers2 == null && rightKmers2 == null) {
                        break;
                    }
                
                }
            }
        }
        
        return new ReadPair(leftKmers, rightKmers, leftCorrected || rightCorrected);
    }
    
    public static String correctErrors(String seq, BloomFilterDeBruijnGraph graph, int lookahead, int errorsAllowed) {
        int numCorrected = 0;
        int seqLen = seq.length();
        int k = graph.getK();
        
        if (seqLen < k) {
            // no correction
            return seq;
        }
        
        int numKmers = seqLen-k+1;
        StringBuilder sb = new StringBuilder(seq);
        
        float bestCov, cov;
        String kmer, guide;
        ArrayDeque<String> variants;
        
        // correct from start
        for (int i=0; i<numKmers; ++i) {
            int j = i+k;
            kmer = sb.substring(i, j);
            variants = graph.getRightVariants(kmer);
            if (!variants.isEmpty()) {
                guide = sb.substring(j, Math.min(j+lookahead, seqLen));
                int guideLen = guide.length();
                bestCov = 0;
                
                if (graph.contains(kmer)) {
                    bestCov = rightGuidedMedianCoverage(graph, kmer, guide, guideLen);
                }
                
                // test for mismatch
                boolean correctedMismatch = false;
                char bestBase = 'N';
                for (String v : variants) {
                    cov = rightGuidedMedianCoverage(graph, v, guide, guideLen);
                                        
                    if (cov > bestCov) {
                        bestCov = cov;
                        bestBase = v.charAt(k-1);
                        correctedMismatch = true;
                    }
                }
                
                // test for insertion in the sequence, ie. last base of kmer is the inserted base
                boolean correctedInsertion = false;
                if (j < seqLen-2) {
                    String a = graph.getPrefix(kmer) + sb.charAt(j);
                    if (graph.contains(a)) {
                        String guideIns = sb.substring(j+1, Math.min(j+1+lookahead, seqLen));

                        cov = rightGuidedMedianCoverage(graph, a, guideIns, guideIns.length());
                        if (cov > bestCov) {
                            correctedMismatch = false;
                            correctedInsertion = true;

                            bestCov = cov;
                        }
                    }
                }
                
                // test for deletion in the sequence
                boolean correctedDeletion = false;
                String guideDel = sb.substring(j-1, Math.min(j-1+lookahead, seqLen));
                int guideDelLen = guideDel.length();
                char bestInsBase = 'N';
                for (String v : variants) {
                    cov = rightGuidedMedianCoverage(graph, v, guideDel, guideDelLen);
                                        
                    if (cov > bestCov) {
                        correctedMismatch = false;
                        correctedInsertion = false;
                        correctedDeletion = true;
                        
                        bestCov = cov;
                        bestInsBase = v.charAt(k-1);
                    }
                }
                

                if (correctedMismatch) {
                    if (++numCorrected > errorsAllowed) {
                        // too many errors; do not apply corrections
                        return seq;
                    }

                    // replace the mismatch base
                    sb.setCharAt(j-1, bestBase);
                }
                else if (correctedInsertion) {
                    if (++numCorrected > errorsAllowed) {
                        // too many errors; do not apply corrections
                        return seq;
                    }

                    // remove the inserted base
                    sb.deleteCharAt(j-1);
                    --seqLen;
                    --numKmers;
                }
                else if (correctedDeletion) {
                    if (++numCorrected > errorsAllowed) {
                        // too many errors; do not apply corrections
                        return seq;
                    }

                    // insert the deleted base
                    sb.insert(j-1, bestInsBase);
                    ++seqLen;
                    ++numKmers;
                }
            }
        }
        
        // correct from end
        for (int i=seqLen-k; i>=0; --i) {
            kmer = sb.substring(i, i+k);
            variants = graph.getLeftVariants(kmer);
            if (!variants.isEmpty()) {
                guide = sb.substring(Math.max(0, i-lookahead), i);
                int guideLen = guide.length();
                bestCov = 0;
                
                if (graph.contains(kmer)) {
                    bestCov = leftGuidedMedianCoverage(graph, kmer, guide, guideLen);
                }
                                
                // test for mismatch
                boolean correctedMismatch = false;
                char bestBase = 'N';
                for (String v : variants) {
                    cov = leftGuidedMedianCoverage(graph, v, guide, guideLen);
                                        
                    if (cov > bestCov) {
                        bestCov = cov;
                        bestBase = v.charAt(0);
                        correctedMismatch = true;
                    }
                }
                
                // test for insertion in the sequence; ie. first base of kmer is the inserted base
                boolean correctedInsertion = false;
                if (i > 0) {
                    String a = sb.charAt(i-1) + graph.getSuffix(kmer);
                    if (graph.contains(a)) {
                        String guideIns = sb.substring(Math.max(0, i-1-lookahead), i-1);

                        cov = leftGuidedMedianCoverage(graph, a, guideIns, guideIns.length());
                        if (cov > bestCov) {
                            correctedMismatch = false;
                            correctedInsertion = true;

                            bestCov = cov;
                        }
                    }
                }
                
                // test for deletion in the sequence
                boolean correctedDeletion = false;
                String guideDel = sb.substring(Math.max(0, i+1-lookahead), i+1);
                int guideDelLen = guideDel.length();
                char bestInsBase = 'N';
                for (String v : variants) {
                    cov = leftGuidedMedianCoverage(graph, v, guideDel, guideDelLen);
                                        
                    if (cov > bestCov) {
                        correctedMismatch = false;
                        correctedInsertion = false;
                        correctedDeletion = true;
                        
                        bestCov = cov;
                        bestInsBase = v.charAt(0);
                    }
                }
                
                if (correctedMismatch) {
                    if (++numCorrected > errorsAllowed) {
                        // too many errors; do not apply corrections
                        return seq;
                    }

                    // replace the mismatch base
                    sb.setCharAt(i, bestBase);
                }
                else if (correctedInsertion) {
                    if (++numCorrected > errorsAllowed) {
                        // too many errors; do not apply corrections
                        return seq;
                    }

                    // remove the inserted base
                    sb.deleteCharAt(i);
                    --seqLen;
                    --numKmers;
                }
                else if (correctedDeletion) {
                    if (++numCorrected > errorsAllowed) {
                        // too many errors; do not apply corrections
                        return seq;
                    }

                    // insert the deleted base
                    sb.insert(i, bestInsBase);
                    ++seqLen;
                    ++numKmers;
                }
            }
        }
        
        if (numCorrected == 0) {
            return seq;
        }
        
        String seq2 = sb.toString();
        if (graph.isValidSeq(seq2)) {
            return seq2;
        }
        
        return seq;
    }
    
    public static float[] coverageGradients(String seq, BloomFilterDeBruijnGraph graph, int lookahead) {
        float[] counts = graph.getCounts(seq);
        int numCounts = counts.length;
        
        ArrayDeque<Float> window = new ArrayDeque<>();
        for (int i=0; i<lookahead; ++i) {
            window.addLast(counts[i]);
        }        

        int numMedCounts = numCounts-lookahead+1;
        float[] medCounts = new float[numMedCounts];
        medCounts[0] = getMedian(window);
        int m = 0;
        for (int i=lookahead; i<numCounts; ++i) {
            window.removeFirst();
            window.addLast(counts[i]);
            medCounts[++m] = getMedian(window);
        }
        
        int numGradients = numCounts-(2*lookahead)+1;
        float[] gradients = new float[numGradients];
        for (int i=0; i<numGradients; ++i) {
            float r = medCounts[i]/medCounts[i+lookahead];
            if (r > 1) {
                gradients[i] = 1/r;
            }
            else {
                gradients[i] = r;
            }
        }
        
        return gradients;
    }

    public static String assemble(ArrayDeque<Kmer> kmers, int k) {
        StringBuilder sb = new StringBuilder(kmers.size() + k - 1);

        int lastIndex = k - 1;
        
        byte[] bytes = kmers.getFirst().bytes;
        for (int i=0; i<lastIndex; ++i) {
            sb.append((char) bytes[i]);
        }
        
        for (Kmer e : kmers) {
            sb.append((char) e.bytes[lastIndex]);
        }
        
        return sb.toString();
    }
    
    public static String assemble(ArrayList<Kmer> kmers, int k, int start, int end) {
        StringBuilder sb = new StringBuilder(end - start + k - 1);

        int lastIndex = k - 1;
        byte[] bytes = kmers.get(start).bytes;
        for (int i=0; i<k; ++i) {
            sb.append((char) bytes[i]);
        }
        
        for (int i=start+1; i<end; ++i) {
            sb.append((char) kmers.get(i).bytes[lastIndex]);
        }
        
        return sb.toString();
    }
    
    public static String assembleReverseOrder(ArrayDeque<Kmer> kmers, int k) {
        StringBuilder sb = new StringBuilder(kmers.size() + k - 1);

        int lastIndex = k - 1;  
        
        Iterator<Kmer> itr = kmers.descendingIterator();
        byte[] bytes = itr.next().bytes;
        for (int i=0; i<k; ++i) {
            sb.append((char) bytes[i]);
        }
        
        while (itr.hasNext()) {
            sb.append((char) itr.next().bytes[lastIndex]);
        }
        
        return sb.toString();
    }
    
    public static String assemble(Collection<Kmer> kmers, int k) {
        StringBuilder sb = new StringBuilder(kmers.size() + k - 1);

        int lastIndex = k - 1;
        
        Iterator<Kmer> itr = kmers.iterator();
        byte[] bytes = itr.next().bytes;
        for (int i=0; i<k; ++i) {
            sb.append((char) bytes[i]);
        }
        
        while (itr.hasNext()) {
            sb.append((char) itr.next().bytes[lastIndex]);
        }

        // surprisingly slower!
//        int numKmers = kmers.size();
//        int length = numKmers + k - 1;
//        StringBuilder sb = new StringBuilder(length);
//
//        for (int i=0; i<numKmers; i+k) {
//            sb.append(kmers.get(i).seq);
//        }
//
//        int remainder = length % k;
//        if (remainder > 0) {
//            sb.append(kmers.get(numKmers-1).seq.substring(k-remainder));
//        }
                
        return sb.toString();
    }

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
    
    public static String assembleFirstBase(Collection<Kmer> kmers, int k) {
        StringBuilder sb = new StringBuilder(kmers.size() + k - 1);
        for (Kmer kmer : kmers) {
            sb.append((char) kmer.bytes[0]);
        }
        
        return sb.toString();
    }

    public static String assembleLastBase(Collection<Kmer> kmers, int k) {
        int lastIndex = k-1;
        
        StringBuilder sb = new StringBuilder(kmers.size() + k - 1);
        for (Kmer kmer : kmers) {
            sb.append((char) kmer.bytes[lastIndex]);
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
        
    public static ArrayList<Kmer> getKmers(String seq, BloomFilterDeBruijnGraph graph, int maxIndelSize, int lookahead) {
        int k = graph.getK();
        int seqLength = seq.length();
        int numKmers = seqLength - k + 1;
        
        ArrayList<Kmer> result = new ArrayList<>(numKmers);
        ArrayList<Kmer> bestResult = null;
        
        byte[] bytes = stringToBytes(seq, seqLength);
        
        NTHashIterator itr = graph.getHashIterator();
        itr.start(seq);
        long[] hVals = itr.hVals;
        int i;
        
        float c;
        int numMissingKmers = 0;
        Kmer nextKmer;
        while (itr.hasNext()) {
            itr.next();
            i = itr.getPos();
            c = graph.getCount(hVals);
            
            if (c > 0) {
                nextKmer = new Kmer(Arrays.copyOfRange(bytes, i, i+k), c, hVals);
                
                if (numMissingKmers > 0 && !result.isEmpty()) {
                    ArrayDeque<Kmer> path = getMaxCoveragePath(graph, result.get(result.size()-1), nextKmer, maxIndelSize + numMissingKmers, lookahead);
                    
                    if (path == null) {
                        if (bestResult == null) {
                            bestResult = result;
                        }
                        else {
                            if (result.size() > bestResult.size()) {
                                bestResult = result;
                            }
                        }
                        
                        result = new ArrayList<>(numKmers);
                    }
                    else {
                        result.addAll(path);
                        result.add(nextKmer);
                        numMissingKmers = 0;
                    }
                }
                else {
                    result.add(nextKmer);
                    numMissingKmers = 0;
                }
            }
            else {
                ++numMissingKmers;
            }
        }
        
        if (bestResult == null || result.size() > bestResult.size()) {
            return result;
        }
        
        return bestResult;
    }
    
    public static ArrayList<Kmer> connectKmers(ArrayList<String> segments, BloomFilterDeBruijnGraph graph, int bound, int lookahead, int maxIndelSize) {
        int numSeqs = segments.size();
        switch (numSeqs) {
            case 0:
                return null;
            case 1:
                return getKmers(segments.get(0), graph, maxIndelSize, lookahead);
            default:
                ArrayList<Kmer> current = getKmers(segments.get(0), graph, maxIndelSize, lookahead);
                ArrayList<Kmer> best = null;
                
                ArrayList<Kmer> next;
                for (int i=1; i<numSeqs; ++i) {
                    next = getKmers(segments.get(i), graph, maxIndelSize, lookahead);
                    
                    if (!next.isEmpty()) {                    
                        ArrayDeque<Kmer> path = getMaxCoveragePath(graph, current.get(current.size()-1), next.get(0), bound, lookahead);

                        if (path == null) {
                            if (best == null) {
                                best = current;
                            }
                            else {
                                if (current.size() > best.size()) {
                                    best = current;
                                }
                            }
                            
                            current = next;
                        }
                        else {
                            current.addAll(path);
                            current.addAll(next);
                        }
                    }
                }
                
                if (best == null || current.size() > best.size()) {
                    return current;
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
        
        ArrayDeque<Kmer> pathKmers = getMaxCoveragePath(graph, graph.getKmer(getLastKmer(left, k)), graph.getKmer(getFirstKmer(right, k)), bound, lookahead);
        
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
        
        return leftWing + assemble(pathKmers, k) + rightWing;
    }
    
    public static ArrayList<Kmer> overlap(ArrayList<Kmer> leftKmers, ArrayList<Kmer> rightKmers, BloomFilterDeBruijnGraph graph, int minOverlap) {
        int k = graph.getK();
        
        return overlap(assemble(leftKmers, k), assemble(rightKmers, k), graph, minOverlap);
    }
    
    public static ArrayList<Kmer> overlap(String left, String right, BloomFilterDeBruijnGraph graph, int minOverlap) {
//        int k = graph.getK();
        
        String overlapped = overlapMaximally(left, right, minOverlap);
        
        if (overlapped != null) {
            ArrayList<Kmer> overlappedKmers = graph.getKmers(overlapped);
            for (Kmer kmer : overlappedKmers) {
                if (kmer.count <= 0) {
                    return null;
                }
            }

            return overlappedKmers;
        }
        
        return null;
    }
    
    public static ArrayList<Kmer> overlapAndConnect(ArrayList<Kmer> leftKmers, ArrayList<Kmer> rightKmers, BloomFilterDeBruijnGraph graph, int bound, int lookahead, int minOverlap) {
        int k = graph.getK();
        String leftSeq = assemble(leftKmers, k);
        String rightSeq = assemble(rightKmers, k);

        // 1. Attempt simple overlap
        ArrayList<Kmer> fragmentKmers = overlap(leftSeq, rightSeq, graph, minOverlap);
        
        if (fragmentKmers == null) {
            Kmer leftLastKmer = leftKmers.get(leftKmers.size()-1);
            Kmer rightFirstKmer = rightKmers.get(0);

            // 2. Attempt connect a path
            ArrayDeque<Kmer> connectedPath = getMaxCoveragePath(graph, leftLastKmer, rightFirstKmer, bound, lookahead);

            if (connectedPath != null) {
                fragmentKmers = new ArrayList<>(leftKmers.size() + connectedPath.size() + rightKmers.size());
                fragmentKmers.addAll(leftKmers);
                fragmentKmers.addAll(connectedPath);
                fragmentKmers.addAll(rightKmers);
            }
//            else {
//                // 3. Attempt dovetail overlap (ie. when fragment is shorter than read length)
//                fragmentKmers = overlap(rightSeq, leftSeq, graph, minOverlap);
//            }
        }
        
        return fragmentKmers;
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
    
    private static LinkedList<Kmer> getSuccessorsRanked(Kmer source, BloomFilterDeBruijnGraph graph, int lookahead) {
        LinkedList<Kmer> results = new LinkedList<>();
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : graph.getSuccessors(source)) {
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
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : graph.getPredecessors(source)) {
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
        final int numKmers = kmers.size();
        
        float minCov = getMinimumKmerCoverage(kmers, Math.max(0, numKmers-pairedKmerDistance), numKmers-1);
        int k = graph.getK();
        int minAnchorDistanceFromEdge = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(minCov), pairedKmerDistance-k);
        
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
            if (!assembledKmersBloomFilter.lookup(kmers.get(start).hashVals)) {
                break;
            }
        }
                
        Kmer kmer, p;
        end += minNumPairs;
        for (int i=Math.min(start, numKmers-minAnchorDistanceFromEdge); i>=end; --i) {
            kmer = kmers.get(i);
            
            if (graph.lookupLeftKmer(kmer.hashVals) &&
                    !isLowComplexity2(kmer.bytes)) {
                int numPaired = 1;
                for (int j=1; j<minNumPairs; ++j) {
                    p = kmers.get(i-j);
                    if (!graph.lookupLeftKmer(p.hashVals) ||
                            isLowComplexity2(p.bytes)) {
                        i=j;
                        break;
                    }
                    else if (++numPaired >= minNumPairs) {
                        return pairedKmerDistance - (numKmers - 1 - i);
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
        
        final int numKmers = kmers.size();
        
        float minCov = getMinimumKmerCoverage(kmers, Math.max(0, numKmers-pairedKmerDistance), numKmers-1);
        int k = graph.getK();
        int minAnchorDistanceFromEdge = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(minCov), pairedKmerDistance-k);
        
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
            if (!assembledKmersBloomFilter.lookup(kmers.get(start).hashVals)) {
                break;
            }
        }
                
        Kmer kmer, p;
        end += minNumPairs;
        for (int i=Math.min(start, numKmers-minAnchorDistanceFromEdge); i>=end; --i) {
            kmer = kmers.get(i);
            
            if (graph.lookupRightKmer(kmer.hashVals) &&
                    !isLowComplexity2(kmer.bytes)) {
                int numPaired = 1;
                for (int j=1; j<minNumPairs; ++j) {
                    p = kmers.get(i-j);
                    if (!graph.lookupRightKmer(p.hashVals) ||
                            isLowComplexity2(p.bytes)) {
                        i=j;
                        break;
                    }
                    else if (++numPaired >= minNumPairs) {
                        return pairedKmerDistance - (numKmers - 1 - i);
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
    
    private static int maxRightPartnerSearchDepth(ArrayList<Kmer> fragmentKmers, BloomFilterDeBruijnGraph graph, int pairedKmerDistance, BloomFilter assembledKmersBloomFilter) {
        
        final int numKmers = fragmentKmers.size();
        
        float[] covs = new float[Math.min(pairedKmerDistance, numKmers)];
        for (int i=0; i<covs.length; ++i) {
            covs[i] = fragmentKmers.get(numKmers-i-1).count;
        }
        
        int k = graph.getK();
        int anchorLength = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(getMinium(covs)), pairedKmerDistance-k);
        if (anchorLength == 0) {
            ++anchorLength;
        }
        
        
        int end = Math.max(0, numKmers-pairedKmerDistance);
        int start;

        for (start=numKmers-1; start>end; --start) {
            if (!assembledKmersBloomFilter.lookup(fragmentKmers.get(start).hashVals)) {
                break;
            }
        }
        
        start = Math.min(numKmers-anchorLength, start);
        
        Kmer kmer;
        for (int i=start; i>end; --i) {
            kmer = fragmentKmers.get(i);
            if (graph.lookupLeftKmer(kmer.hashVals) && i>0 && !isLowComplexity2(kmer.bytes)) {
                kmer = fragmentKmers.get(i-1);
                if (graph.lookupLeftKmer(kmer.hashVals) && !isLowComplexity2(kmer.bytes)) {
                    return pairedKmerDistance - (numKmers - i);
                }
                else {
                    --i;
                }
            }
        }
        
        return 0;
    }
    
    private static int maxLeftPartnerSearchDepth(ArrayList<Kmer> fragmentKmers, BloomFilterDeBruijnGraph graph, int pairedKmerDistance, BloomFilter assembledKmersBloomFilter) {
        
        final int numKmers = fragmentKmers.size();
        
        float[] covs = new float[Math.min(pairedKmerDistance, numKmers)];
        for (int i=0; i<covs.length; ++i) {
            covs[i] = fragmentKmers.get(numKmers-i-1).count;
        }
        
        int k = graph.getK();
        int anchorLength = Math.min(k * RNABloom.getMinCoverageOrderOfMagnitude(getMinium(covs)), pairedKmerDistance-k);
        if (anchorLength == 0) {
            ++anchorLength;
        }
        
        int end = Math.max(0, numKmers-pairedKmerDistance);
        int start;
        
        for (start=numKmers-1; start>end; --start) {
            if (!assembledKmersBloomFilter.lookup(fragmentKmers.get(start).hashVals)) {    
                break;
            }
        }
        
        start = Math.min(numKmers-anchorLength, start);
        
        Kmer kmer;
        for (int i=start; i>end; --i) {
            kmer = fragmentKmers.get(i);
            if (graph.lookupRightKmer(kmer.hashVals) && i>0 && !isLowComplexity2(kmer.bytes)) {
                kmer = fragmentKmers.get(i-1);
                if (graph.lookupRightKmer(kmer.hashVals) && !isLowComplexity2(kmer.bytes)) {
                    return pairedKmerDistance - (numKmers - i);
                }
                else {
                    --i;
                }
            }
        }
        
        return 0;
    }
    
    private static boolean hasPairedRightKmers(Kmer source,
                                            ArrayList<Kmer> kmers,
                                            int partnerFromIndex,
                                            int partnerToIndex,
                                            BloomFilterDeBruijnGraph graph) {
        Kmer kmer, partner;
        
        partner = kmers.get(partnerFromIndex);
        
        if (graph.lookupKmerPair(partner.hashVals, source.hashVals)) {
            
            ArrayDeque<Kmer> frontier = graph.getSuccessors(source);
            ArrayDeque<Kmer> newFrontier = new ArrayDeque();
            ArrayDeque<Kmer> tmp;

            Iterator<Kmer> itr;

            for (int i=partnerFromIndex+1; i<partnerToIndex; ++i) {
                partner = kmers.get(i);

                if (!graph.lookupLeftKmer(partner.hashVals)) {
                    return false;
                }

                itr = frontier.iterator();
                while (itr.hasNext()) {
                    kmer = itr.next();
                    if (graph.lookupRightKmer(kmer.hashVals) && 
                            graph.lookupKmerPairing(partner.hashVals, kmer.hashVals)) {
                        newFrontier.addAll(graph.getSuccessors(kmer));
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
        
        Kmer kmer, partner;
        
        partner = kmers.get(partnerFromIndex);
        
        if (graph.lookupKmerPair(source.hashVals, partner.hashVals)) {
            
            ArrayDeque<Kmer> frontier = graph.getPredecessors(source);
            ArrayDeque<Kmer> newFrontier = new ArrayDeque<>();
            ArrayDeque<Kmer> tmp;

            Iterator<Kmer> itr;

            for (int i=partnerFromIndex+1; i<partnerToIndex; ++i) {
                partner = kmers.get(i);

                if (!graph.lookupRightKmer(partner.hashVals)) {
                    return false;
                }

                itr = frontier.iterator();
                while (itr.hasNext()) {
                    kmer = itr.next();
                    if (graph.lookupLeftKmer(kmer.hashVals) && 
                            graph.lookupKmerPairing(kmer.hashVals, partner.hashVals)) {
                        newFrontier.addAll(graph.getPredecessors(kmer));
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
            if (!graph.lookupLeftKmer(kmers.get(i).hashVals)) {
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
            if (!graph.lookupRightKmer(kmers.get(i).hashVals)) {
                return false;
            }
        }
        
        return true;
    }
    
    private static boolean extendRightWithPairedKmersBFS(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs,
                                            BloomFilter assembledKmersBloomFilter,
                                            HashSet<String> usedKmers,
                                            HashSet<String> usedPartnerKmers) {
        
        int distance = graph.getPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int numKmers = kmers.size();
                            
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = graph.getSuccessors(cursor);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        Iterator<Kmer> itr;
        int partnerIndex = numKmers-distance;
        Kmer kmer, partner;
        ArrayDeque<Kmer> simpleExtension;
        
        while (!neighbors.isEmpty()) {
            simpleExtension = extendRight(cursor,
                                            graph, 
                                            maxTipLength, 
                                            usedKmers, 
                                            maxIndelSize, 
                                            percentIdentity);
                        
            if (!simpleExtension.isEmpty()) {
                itr = simpleExtension.descendingIterator();
                while(itr.hasNext() && assembledKmersBloomFilter.lookup(itr.next().hashVals)) {
                    itr.remove();
                }
                
                if (!simpleExtension.isEmpty()) {
                    kmers.addAll(simpleExtension);
                    numKmers = kmers.size();

                    if (numKmers < distance) {
                        // too short to have partners
                        return true;
                    }

                    partnerIndex += simpleExtension.size();
                    
                    for (Kmer e : simpleExtension) {
                        usedKmers.add(e.toString());
                    }

                    // NOTE: kmer at `partnerIndex` will be paired with `cursor`
                    for (int i=Math.max(0, partnerIndex-simpleExtension.size()); i<partnerIndex; ++i) {
                        usedPartnerKmers.add(kmers.get(i).toString());
                    }
                    
                    neighbors = graph.getSuccessors(simpleExtension.getLast());
                    if (neighbors.isEmpty()) {
                        return false;
                    }
                }
            }
            
            if (numKmers < distance) {
                // too short to have partners
                return true;
            }
            
            cursor = null;
                        
            if (!areLeftKmers(kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                // not enough supporting partners
                return true;
            }
            
            float minEdgeCoverage = getMinimumKmerCoverage(kmers, numKmers-distance, numKmers-1);
            float minCovThreshold = minEdgeCoverage * 0.5f;
            ArrayDeque<Kmer> pairedNeighbors = new ArrayDeque<>(4);
            itr = neighbors.iterator();
            while (itr.hasNext()) {
                kmer = itr.next();
                
                if (hasPairedRightKmers(kmer, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                    pairedNeighbors.add(kmer);
                    
                    if (kmer.count < minCovThreshold) {
                        itr.remove();
                    }
                }
                else {
                    itr.remove();
                }
            }
            
            if (neighbors.isEmpty()) {
                Kmer highestCoverageNeighbor = greedyExtendRightOnce(graph, pairedNeighbors, lookahead);
                if (highestCoverageNeighbor == null) {
                    return true;
                }
                else {
                    cursor = highestCoverageNeighbor;
                }
            }
            else if (neighbors.size() == 1) {
                cursor = neighbors.pop();
            }
            else {
                float minDiffCoverage = Float.MAX_VALUE;
                
                float c,d;
                for (Kmer n : neighbors) {
                    c = getMaxMedianCoverageRight(graph, n, lookahead);
                    if (c > 0) {
                        d = Math.abs(minEdgeCoverage - c);
                        if (d < minDiffCoverage) {
                            minDiffCoverage = d;
                            cursor = n;
                        }
                    }
                }
                
//                cursor = greedyExtendRightOnce(graph, neighbors, lookahead);
//                
//                if (cursor == null) {
//                    // no good candidates
//                    return true;
//                }
            }
            
            partner = kmers.get(partnerIndex);
            
            String cursorSeq = cursor.toString();
            String partnerSeq = partner.toString();
            
            if (usedKmers.contains(cursorSeq) &&
                    usedPartnerKmers.contains(partnerSeq)){
                // loop, do not extend further
                return false;
            }
            
            if (assembledKmersBloomFilter.lookup(cursor.hashVals)) {
                boolean assembled = greedyExtendRight(graph, cursor, lookahead, lookahead, assembledKmersBloomFilter) != null;
                
                if (assembled) {
//                    int numNotAssembled = 0;
                    for (int i=partnerIndex; i<numKmers; ++i) {
                        if (!assembledKmersBloomFilter.lookup(kmers.get(i).hashVals)) {
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
            
            usedKmers.add(cursorSeq);
            usedPartnerKmers.add(partnerSeq);
            
            ++partnerIndex;
            ++numKmers;
            
            neighbors = graph.getSuccessors(cursor);
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
                                            HashSet<String> usedKmers,
                                            HashSet<String> usedPartnerKmers) {
        
        int distance = graph.getPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int numKmers = kmers.size();
        
        // Note that `kmers` are in reverse order already
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = graph.getPredecessors(cursor);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        Iterator<Kmer> itr;
        int partnerIndex = numKmers-distance;
        Kmer kmer, partner;
        ArrayDeque<Kmer> simpleExtension;
        
        while (!neighbors.isEmpty()) {
            simpleExtension = extendLeft(cursor,
                                            graph, 
                                            maxTipLength, 
                                            usedKmers, 
                                            maxIndelSize, 
                                            percentIdentity);
            
            if (!simpleExtension.isEmpty()) {
                itr = simpleExtension.descendingIterator();
                while(itr.hasNext() && assembledKmersBloomFilter.lookup(itr.next().hashVals)) {
                    itr.remove();
                }
                
                if (!simpleExtension.isEmpty()) {
                    kmers.addAll(simpleExtension);
                    numKmers = kmers.size();

                    if (numKmers < distance) {
                        // too short to have partners
                        return true;
                    }

                    partnerIndex += simpleExtension.size();

                    for (Kmer e : simpleExtension) {
                        usedKmers.add(e.toString());
                    }

                    // NOTE: kmer at `partnerIndex` will be paired with `cursor`
                    for (int i=Math.max(0, partnerIndex-simpleExtension.size()); i<partnerIndex; ++i) {
                        usedPartnerKmers.add(kmers.get(i).toString());
                    }
                    
                    neighbors = graph.getPredecessors(simpleExtension.getLast());
                    if (neighbors.isEmpty()) {
                        return false;
                    }
                }
            }
            
            if (numKmers < distance) {
                // too short to have partners
                return true;
            }
            
            cursor = null;
            
            if (!areRightKmers(kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                // not enough supporting partners
                return true;
            }
            
            float minEdgeCoverage = getMinimumKmerCoverage(kmers, numKmers-distance, numKmers-1);
            float minCovThreshold = minEdgeCoverage * 0.5f;
            ArrayDeque<Kmer> pairedNeighbors = new ArrayDeque<>(4);
            itr = neighbors.iterator();
            while (itr.hasNext()) {
                kmer = itr.next();
                
                if (hasPairedLeftKmers(kmer, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                    pairedNeighbors.add(kmer);
                    
                    if (kmer.count < minCovThreshold) {
                        itr.remove();
                    }
                }
                else {
                    itr.remove();
                }
            }
            
            if (neighbors.isEmpty()) {
                Kmer highestCoverageNeighbor = greedyExtendLeftOnce(graph, pairedNeighbors, lookahead);
                if (highestCoverageNeighbor == null) {
                    return true;
                }
                else {
                    cursor = highestCoverageNeighbor;
                }
            }
            else if (neighbors.size() == 1) {
                cursor = neighbors.pop();
            }
            else {
                float minDiffCoverage = Float.MAX_VALUE;
                
                float c,d;
                for (Kmer n : neighbors) {
                    c = getMaxMedianCoverageLeft(graph, n, lookahead);
                    if (c > 0) {
                        d = Math.abs(minEdgeCoverage - c);
                        if (d < minDiffCoverage) {
                            minDiffCoverage = d;
                            cursor = n;
                        }
                    }
                }
                
//                cursor = greedyExtendLeftOnce(graph, neighbors, lookahead);
//                
//                if (cursor == null) {
//                    // no good candidates
//                    return true;
//                }
            }
            
            partner = kmers.get(partnerIndex);
            
            String cursorSeq = cursor.toString();
            String partnerSeq = partner.toString();
            
            if (usedKmers.contains(cursorSeq) &&
                    usedPartnerKmers.contains(partnerSeq)){
                // loop, do not extend further
                return false;
            }
            
            if (assembledKmersBloomFilter.lookup(cursor.hashVals)) {
                boolean assembled = greedyExtendLeft(graph, cursor, lookahead, lookahead, assembledKmersBloomFilter) != null;
                
                if (assembled) {
//                    int numNotAssembled = 0;
                    for (int i=partnerIndex; i<numKmers; ++i) {
                        if (!assembledKmersBloomFilter.lookup(kmers.get(i).hashVals)) {
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
            
            usedKmers.add(cursorSeq);
            usedPartnerKmers.add(partnerSeq);
                        
            ++partnerIndex;
            ++numKmers;
            
            neighbors = graph.getPredecessors(cursor);
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
                                            HashSet<String> usedKmers,
                                            HashSet<String> usedPartnerKmers) {
        
        final int distance = graph.getPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int maxDepth = maxRightPartnerSearchDepth2(kmers, graph, distance, assembledKmersBloomFilter, minNumPairs);
        
        // data structure to store visited kmers at defined depth
        HashMap<String, ArrayDeque<Integer>> visitedKmers = new HashMap();
                
        int numKmers = kmers.size();
        int depth = 0;
        int partnerIndex = numKmers - distance + depth;
        
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
        
        branchesStack.add(getSuccessorsRanked(kmers.get(numKmers-1), graph, lookahead));
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        HashSet<String> extensionKmers = new HashSet<>();
        
        Kmer cursor;
        int maxPartnerIndex = numKmers - 1 - minNumPairs;
        while (!branchesStack.isEmpty()) {
            LinkedList<Kmer> branches = branchesStack.getLast();
            
            if (branches.isEmpty()) {
                cursor = extension.pollLast();
                
                if (cursor != null) {
                    extensionKmers.remove(cursor.toString());
                }
                
                branchesStack.removeLast();
                --depth;
                --partnerIndex;
            }
            else {
                cursor = branches.pop();
                String cursorSeq = cursor.toString();
                
                if (partnerIndex >=0 &&
                        partnerIndex <= maxPartnerIndex &&
//                        partnerIndex+minNumPairs < kmers.size() && 
                        hasPairedRightKmers(cursor, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                    
                    String partnerSeq = kmers.get(partnerIndex).toString();
                    
                    if (usedKmers.contains(cursorSeq) && usedPartnerKmers.contains(partnerSeq)) {
                        return false;
                    }
                                        
                    if (assembledKmersBloomFilter.lookup(cursor.hashVals)) {
                        boolean assembled = greedyExtendRight(graph, cursor, lookahead, lookahead, assembledKmersBloomFilter) != null;

                        if (assembled) {
//                            int numNotAssembled = 0;
                            
                            for (Kmer kmer : extension) {
                                if (!assembledKmersBloomFilter.lookup(kmer.hashVals)) {
//                                    if (distanceInversePI < ++numNotAssembled) {
                                        assembled = false;
                                        break;
//                                    }
                                }
                            }
                            
                            if (assembled) {
                                for (int i=partnerIndex; i<numKmers; ++i) {
                                    if (!assembledKmersBloomFilter.lookup(kmers.get(i).hashVals)) {
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
                    
                    for (Kmer e : extension) {
                        usedKmers.add(e.toString());
                    }
                    usedKmers.add(cursorSeq);
                                        
                    for (int i=Math.max(0, partnerIndex-extension.size()); i<=partnerIndex; ++i) {
                        usedPartnerKmers.add(kmers.get(i).toString());
                    }
                    
                    return true;
                }
                else if (depth < maxDepth &&
                        (depth == 0 || !extensionKmers.contains(cursorSeq))) {
                    
                    if (graph.hasAtLeastXPredecessors(cursor, 2)) {
                        // only consider kmers that may be visited from an alternative branch upstream
                        
                        ArrayDeque<Integer> visitedDepths = visitedKmers.get(cursorSeq);
                        if (visitedDepths == null) {
                            visitedDepths = new ArrayDeque<>();
                            visitedDepths.add(depth);
                            
                            visitedKmers.put(cursorSeq, visitedDepths);
                            
                            branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursorSeq);
                            ++depth;
                            ++partnerIndex;
                        }
                        else {
                            boolean visited = false;
                            for (Integer d : visitedDepths) {
                                if (d == depth) {
                                    visited = true;
                                    break;
                                }
                            }
                            
                            if (!visited) {
                                visitedDepths.add(depth);

                                branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursorSeq);
                                ++depth;
                                ++partnerIndex;
                            }
                        }
                    }
                    else {
                        branchesStack.add(getSuccessorsRanked(cursor, graph, lookahead));
                        extension.add(cursor);
                        extensionKmers.add(cursorSeq);
                        ++depth;
                        ++partnerIndex;
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
                                            HashSet<String> usedKmers,
                                            HashSet<String> usedPartnerKmers) {
        
        final int distance = graph.getPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int maxDepth = maxLeftPartnerSearchDepth2(kmers, graph, distance, assembledKmersBloomFilter, minNumPairs);
        
        // data structure to store visited kmers at defined depth
        HashMap<String, ArrayDeque<Integer>> visitedKmers = new HashMap();
        
        int numKmers = kmers.size();
        int depth = 0;
        int partnerIndex = numKmers - distance + depth;
        
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
        
        branchesStack.add(getPredecessorsRanked(kmers.get(numKmers-1), graph, lookahead));
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        HashSet<String> extensionKmers = new HashSet<>();
        
        Kmer cursor;
        int maxPartnerIndex = numKmers - 1 - minNumPairs;
        while (!branchesStack.isEmpty()) {
            LinkedList<Kmer> branches = branchesStack.getLast();
            
            if (branches.isEmpty()) {
                cursor = extension.pollLast();
                if (cursor != null) {
                    extensionKmers.remove(cursor.toString());
                }
                
                branchesStack.removeLast();
                --depth;
                --partnerIndex;
            }
            else {
                cursor = branches.pop();
                String cursorSeq = cursor.toString();
                
                if (partnerIndex >=0 &&
                        partnerIndex <= maxPartnerIndex &&
//                        partnerIndex+minNumPairs < kmers.size() && 
                        hasPairedLeftKmers(cursor, kmers, partnerIndex, partnerIndex+minNumPairs, graph)) {
                    
                    String partnerSeq = kmers.get(partnerIndex).toString();
                    
                    if (usedKmers.contains(cursorSeq) && usedPartnerKmers.contains(partnerSeq)) {
                        return false;
                    }
                                        
                    if (assembledKmersBloomFilter.lookup(cursor.hashVals)) {
                        boolean assembled = greedyExtendLeft(graph, cursor, lookahead, lookahead, assembledKmersBloomFilter) != null;

                        if (assembled) {
//                            int numNotAssembled = 0;
                            
                            for (Kmer kmer : extension) {
                                if (!assembledKmersBloomFilter.lookup(kmer.hashVals)) {
//                                    if (distanceInversePI < ++numNotAssembled) {
                                        assembled = false;
                                        break;
//                                    }
                                }
                            }
                            
                            if (assembled) {
                                for (int i=partnerIndex; i<numKmers; ++i) {
                                    if (!assembledKmersBloomFilter.lookup(kmers.get(i).hashVals)) {
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
                    
                    for (Kmer e : extension) {
                        usedKmers.add(e.toString());
                    }
                    usedKmers.add(cursorSeq);
                    
                    for (int i=Math.max(0, partnerIndex-extension.size()); i<=partnerIndex; ++i) {
                        usedPartnerKmers.add(kmers.get(i).toString());
                    }
                    
                    return true;
                }
                else if (depth < maxDepth &&
                        (depth == 0 || !extensionKmers.contains(cursorSeq))) {
                    
                    if (graph.hasAtLeastXSuccessors(cursor, 2)) {
                        // only consider kmers that may be visited from an alternative branch upstream
                        
                        ArrayDeque<Integer> visitedDepths = visitedKmers.get(cursorSeq);
                        if (visitedDepths == null) {
                            visitedDepths = new ArrayDeque<>();
                            visitedDepths.add(depth);
                            
                            visitedKmers.put(cursorSeq, visitedDepths);
                            
                            branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                            extension.add(cursor);
                            extensionKmers.add(cursorSeq);
                            ++depth;
                            ++partnerIndex;
                        }
                        else {
                            boolean visited = false;
                            for (Integer d : visitedDepths) {
                                if (d == depth) {
                                    visited = true;
                                    break;
                                }
                            }
                            
                            if (!visited) {
                                visitedDepths.add(depth);

                                branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                                extension.add(cursor);
                                extensionKmers.add(cursorSeq);
                                ++depth;
                                ++partnerIndex;
                            }
                        }
                    }
                    else {
                        branchesStack.add(getPredecessorsRanked(cursor, graph, lookahead));
                        extension.add(cursor);
                        extensionKmers.add(cursorSeq);
                        ++depth;
                        ++partnerIndex;
                    }
                }
            }
        }
        
        return false;
    }
    
    public static void extendWithPairedKmers2(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength,
                                            BloomFilter assembledKmersBloomFilter,
                                            int maxIndelSize,
                                            float percentIdentity,
                                            int minNumPairs) {
        
        HashSet<String> usedKmers = new HashSet<>();
        for (Kmer kmer : kmers) {
            usedKmers.add(kmer.toString());
        }
        
        // naive extend RIGHT
//        ArrayDeque<Kmer> simpleExtension = extendRight(kmers.get(kmers.size()-1),
//                                                        graph, 
//                                                        maxTipLength, 
//                                                        usedKmers, 
//                                                        maxIndelSize, 
//                                                        percentIdentity);
//        if (!simpleExtension.isEmpty()) {
//            Iterator<Kmer> itr = simpleExtension.descendingIterator();
//            while(itr.hasNext() && assembledKmersBloomFilter.lookup(itr.next().hashVals)) {
//                itr.remove();
//            }
//            
//            kmers.addAll(simpleExtension);
//            for (Kmer e : simpleExtension) {
//                usedKmers.add(e.toString());
//            }
//        }
        
        // naive extend LEFT
        Collections.reverse(kmers);
//        simpleExtension = extendLeft(kmers.get(kmers.size()-1),
//                                        graph, 
//                                        maxTipLength, 
//                                        usedKmers, 
//                                        maxIndelSize, 
//                                        percentIdentity);
//        if (!simpleExtension.isEmpty()) {
//            Iterator<Kmer> itr = simpleExtension.descendingIterator();
//            while(itr.hasNext() && assembledKmersBloomFilter.lookup(itr.next().hashVals)) {
//                itr.remove();
//            }
//            
//            kmers.addAll(simpleExtension);
//            for (Kmer e : simpleExtension) {
//                usedKmers.add(e.toString());
//            }
//        }
        
        int distance = graph.getPairedKmerDistance();
        
        // extend with paired kmers LEFT
        
        HashSet<String> usedPartnerKmers = new HashSet<>();
        for (int i=kmers.size()-1-distance; i>=0; --i) {
            usedPartnerKmers.add(kmers.get(i).toString());
        }
        
        boolean extendable = true;
        
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
                                        usedPartnerKmers);
            
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
                                        usedPartnerKmers);
            }
            else {
                break;
            }
        }
        
        Collections.reverse(kmers);
        
        // extend with paired kmers RIGHT
        
        usedPartnerKmers.clear();
        for (int i=kmers.size()-1-distance; i>=0; --i) {
            usedPartnerKmers.add(kmers.get(i).toString());
        }
        
        extendable = true;
        
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
                                        usedPartnerKmers);
            
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
                                        usedPartnerKmers);
            }
            else {
                break;
            }
        }
    }
   
/*
    public static void extendWithPairedKmers(ArrayList<Kmer> kmers, 
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
//        ArrayList<Kmer> naiveExtension = naiveExtendLeft(kmers.get(0), graph, maxTipLength, usedKmers, true);
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
  
        Kmer kmer, n, n2, partner, partner2;
        int partnerIndex;
        String mergedSeq;
        
        ArrayDeque<Kmer> extension = new ArrayDeque<>();
        
        ArrayDeque<LinkedList<Kmer>> branchesStack = new ArrayDeque<>();
        LinkedList<Kmer> neighbors = getSuccessorsRanked(kmers.get(kmers.size()-1), graph, lookahead);
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

                            ListIterator<Kmer> nItr = neighbors.listIterator();
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
                                    
                                    Iterator<Kmer> kmerItr = extension.iterator();
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

                                    ArrayDeque<Kmer> extendRight = extendRight(n, graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
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
                            
                            Iterator<Kmer> kmerItr = extension.iterator();
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
                            
                            ArrayDeque<Kmer> extendRight = extendRight(n, graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
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
            ArrayDeque<Kmer> naiveExtension = extendRight(kmers.get(kmers.size()-1), graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
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

                            ListIterator<Kmer> nItr = neighbors.listIterator();
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
                                    
                                    Iterator<Kmer> kmerItr = extension.iterator();
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

                                    ArrayDeque<Kmer> extendLeft = extendLeft(n, graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
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
                            
                            Iterator<Kmer> kmerItr = extension.iterator();
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

                            ArrayDeque<Kmer> extendLeft = extendLeft(n, graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
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
            ArrayDeque<Kmer> naiveExtension = extendLeft(kmers.get(kmers.size()-1), graph, maxTipLength, usedKmers, maxIndelSize, percentIdentity);
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
        ArrayDeque<ArrayDeque> frontier = new ArrayDeque<>();
        ArrayDeque<Kmer> alts = graph.getSuccessors(source);
        frontier.add(alts);
        
        while (!frontier.isEmpty()) {
            alts = frontier.peekLast();
            if (alts.isEmpty()) {
                frontier.removeLast();
            }
            else {
                frontier.add(graph.getSuccessors(alts.pop()));
            }

            if (frontier.size() >= depth) {
                return true;
            }
        }
        
        return false;
    }
    
    public static boolean hasDepthLeft(Kmer source, BloomFilterDeBruijnGraph graph, int depth) {
        ArrayDeque<ArrayDeque> frontier = new ArrayDeque<>();
        ArrayDeque<Kmer> alts = graph.getPredecessors(source);
        frontier.add(alts);
        
        while (!frontier.isEmpty()) {
            alts = frontier.peekLast();
            if (alts.isEmpty()) {
                frontier.removeLast();
            }
            else {
                frontier.add(graph.getPredecessors(alts.pop()));
            }

            if (frontier.size() >= depth) {
                return true;
            }
        }
        
        return false;
    }
    
    public static ArrayDeque<Kmer> naiveExtendRight(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<String> terminators) {        
        int k = graph.getK();
        
        HashSet<String> usedKmers = new HashSet<>();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = graph.getSuccessors(kmer);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : graph.getLeftVariants(best)) {
                if (hasDepthLeft(s, graph, maxTipLength)) {
                    return result;
                }
            }
            
            if (neighbors.size() == 1) {
                best = neighbors.peek();
            }
            else {
                best = null;
                for (Kmer n : neighbors) {
                    if (hasDepthRight(n, graph, maxTipLength)) {
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
            
            String bestSeq = best.toString();
            
            if (terminators.contains(bestSeq) || usedKmers.contains(bestSeq)) {
                break;
            }
            
            result.add(best);
            usedKmers.add(bestSeq);
            
            neighbors = graph.getSuccessors(best);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> naiveExtendLeft(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<String> terminators) {        
        int k = graph.getK();
        
        HashSet<String> usedKmers = new HashSet<>();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = graph.getPredecessors(kmer);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : graph.getRightVariants(best)) {
                if (hasDepthRight(s, graph, maxTipLength)) {
                    return result;
                }
            }
            
            if (neighbors.size() == 1) {
                best = neighbors.peek();
            }
            else {
                best = null;
                for (Kmer n : neighbors) {
                    if (hasDepthLeft(n, graph, maxTipLength)) {
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
            
            String bestSeq = best.toString();
            
            if (terminators.contains(bestSeq) || usedKmers.contains(bestSeq)) {
                break;
            }
            
            result.addLast(best);
            usedKmers.add(bestSeq);
            
            neighbors = graph.getPredecessors(best);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> extendRight(Kmer source,
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLength, 
                                            HashSet<String> terminators, 
                                            int maxIndelSize, 
                                            float percentIdentity) {
        
        HashSet<String> usedKmers = new HashSet<>();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = graph.getSuccessors(source);
        Kmer best;
//        HashSet<String> backBranchesVisited = new HashSet<>(3);
        
        int k = graph.getK();
        
        while (!neighbors.isEmpty()) {
//            /** look for back branches*/
//            for (Kmer s : graph.getLeftVariants(best)) {
//                if (!backBranchesVisited.contains(s.toString()) && hasDepthLeft(s, graph, maxTipLength)) {
//                    return result;
//                }
//            }
//            
//            backBranchesVisited.clear();
            
            if (neighbors.size() == 1) {
                best = neighbors.peek();
                
                String bestSeq = best.toString();
                if (terminators.contains(bestSeq) || usedKmers.contains(bestSeq)) {
                    return result;
                }
                else {
                    result.add(best);
                    usedKmers.add(bestSeq);
                }
                
                ArrayDeque<Kmer> b = naiveExtendRight(best, graph, maxTipLength, terminators);
                
                if (!b.isEmpty()) {
                    for (Kmer kmer : b) {
                        String kmerSeq = kmer.toString();
                        if (usedKmers.contains(kmerSeq)) {
                            return result;
                        }
                        
                        result.add(kmer);
                        usedKmers.add(kmerSeq);
                    }
                    best = result.peekLast();
                }
            }
            else {
                best = null;
                
                ArrayDeque<ArrayDeque<Kmer>> branches = new ArrayDeque<>(4);
                
                float maxCov = -1;
                for (Kmer n : neighbors) {
                    if (hasDepthRight(n, graph, maxTipLength)) {
                        ArrayDeque<Kmer> b = naiveExtendRight(n, graph, maxTipLength, terminators);
                        
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
                    
                    for (Kmer n : branches.peekLast()) {
                        String nSeq = n.toString();
                        
                        if (usedKmers.contains(nSeq)) {
                            return result;
                        }
                        
                        result.add(n);
                        usedKmers.add(nSeq);
                    }
                    
                    best = result.peekLast();
                }
                else {
                    // multiple good branches
                    
                    // check whether the branches form a bubble
                    ArrayDeque<Kmer> bestBranch = branches.pollFirst();
                    String bestBranchSeq = assemble(bestBranch, k);
                    String suffix = graph.getSuffix(bestBranch.peekLast().toString());
                    int bestBranchLength = bestBranch.size();
                    
                    for (ArrayDeque<Kmer> b : branches) {
                        int len = b.size();
                        
                        if (len >= bestBranchLength - maxIndelSize && len <= bestBranchLength + maxIndelSize) {
                            // length is within range
                            
                            if (!suffix.equals(graph.getSuffix(b.peekLast().toString())) || 
                                    getPercentIdentity(assemble(b, k), bestBranchSeq) < percentIdentity) {
                                return result;
                            }
                        }
                        else {
                            if (len < bestBranchLength - maxIndelSize) {
                                // compare percent identity
                                if (getPercentIdentity(assemble(b, k), bestBranchSeq.substring(0, len+k-1)) < percentIdentity) {
                                    return result;
                                }
                            }
                            else {
                                return result;
                            }
                        }
                    }
                                        
                    for (Kmer n : bestBranch) {
                        String nSeq = n.toString();
                        if (usedKmers.contains(nSeq)) {
                            return result;
                        }
                        
                        result.add(n);
                        usedKmers.add(nSeq);
                    }
                    
                    best = bestBranch.peekLast();
                    neighbors = graph.getSuccessors(best);
                    if (neighbors.size() == 1) {
                        // bubble branches converge at this kmer
                        best = neighbors.pop();
                        result.add(best);
                        usedKmers.add(best.toString());
                    }
                    
//                    for (ArrayDeque<Kmer> b : branches) {
//                        backBranchesVisited.add(b.peekLast().toString());
//                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            neighbors = graph.getSuccessors(best);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> extendLeft(Kmer source,
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLength, 
                                            HashSet<String> terminators, 
                                            int maxIndelSize, 
                                            float percentIdentity) {
        
        HashSet<String> usedKmers = new HashSet<>();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = graph.getPredecessors(source);
        Kmer best;
//        HashSet<String> backBranchesVisited = new HashSet<>(3);
        
        int k = graph.getK();
        
        while (!neighbors.isEmpty()) {
//            /** look for back branches*/
//            for (Kmer s : graph.getRightVariants(best)) {
//                if (!backBranchesVisited.contains(s.toString()) && hasDepthRight(s, graph, maxTipLength)) {
//                    return result;
//                }
//            }
//            
//            backBranchesVisited.clear();
            
            if (neighbors.size() == 1) {
                best = neighbors.peek();
                String bestSeq = best.toString();
                
                if (terminators.contains(bestSeq) || usedKmers.contains(bestSeq)) {
                    return result;
                }
                
                result.add(best);
                
                ArrayDeque<Kmer> b = naiveExtendLeft(best, graph, maxTipLength, terminators);
                
                if (!b.isEmpty()) {
                    for (Kmer kmer : b) {
                        String kmerSeq = kmer.toString();
                        if (usedKmers.contains(kmerSeq)) {
                            return result;
                        }
                        
                        result.add(kmer);
                        usedKmers.add(kmerSeq);
                    }
                    
                    best = result.peekLast();
                }
            }
            else {
                best = null;
                
                ArrayDeque<ArrayDeque<Kmer>> branches = new ArrayDeque<>(4);
                
                float maxCov = -1;
                for (Kmer n : neighbors) {
                    if (hasDepthLeft(n, graph, maxTipLength)) {
                        ArrayDeque<Kmer> b = naiveExtendLeft(n, graph, maxTipLength, terminators);
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
                    
                    for (Kmer n : branches.peekLast()) {
                        String nSeq = n.toString();
                        if (usedKmers.contains(nSeq)) {
                            return result;
                        }
                        
                        result.add(n);
                        usedKmers.add(nSeq);
                    }
                    
                    best = result.peekLast();
                }
                else {
                    // multiple good branches
                    
                    // check whether the branches form a bubble
                    ArrayDeque<Kmer> bestBranch = branches.pollFirst();
                    String bestBranchSeq = assembleReverseOrder(bestBranch, k);
                    String prefix = graph.getPrefix(bestBranch.peekLast().toString());
                    int bestBranchLength = bestBranch.size();
                    
                    for (ArrayDeque<Kmer> b : branches) {
                        int len = b.size();
                        
                        if (len >= bestBranchLength - maxIndelSize && len <= bestBranchLength + maxIndelSize) {
                            // length is within range
                            
                            if (!prefix.equals(graph.getPrefix(b.peekLast().toString())) || 
                                    getPercentIdentity(assembleReverseOrder(b, k), bestBranchSeq) < percentIdentity) {
                                return result;
                            }
                        }
                        else {
                            if (len < bestBranchLength - maxIndelSize) {
                                // compare percent identity
                                if (getPercentIdentity(assembleReverseOrder(b, k), bestBranchSeq.substring(bestBranchSeq.length() - (len+k-1))) < percentIdentity) {
                                    return result;
                                }
                            }
                            else {
                                return result;
                            }
                        }
                    }
                    
                    for (Kmer n : bestBranch) {
                        String nSeq = n.toString();
                        if (usedKmers.contains(nSeq)) {
                            return result;
                        }
                        
                        result.add(n);
                        usedKmers.add(nSeq);
                    }
                    
                    best = bestBranch.peekLast();
                    neighbors = graph.getPredecessors(best);
                    if (neighbors.size() == 1) {
                        // bubble branches converge at this kmer
                        best = neighbors.pop();
                        result.add(best);
                        usedKmers.add(best.toString());
                    }
                    
//                    for (ArrayDeque<Kmer> b : branches) {
//                        backBranchesVisited.add(b.peekLast().toString());
//                    }
                }
            }
            
            if (best == null) {
                break;
            }
            
            neighbors = graph.getPredecessors(best);
        }
        
        return result;
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
//        ArrayList<Kmer> kmers = new ArrayList<>();
//        for (String kmer : SeqUtils.kmerize(seq, k)) {
//             kmers.add(new Kmer(kmer, 1, new long[0]));
//        }
//
//        System.out.println(assemble(kmers, k));
//    }
}
