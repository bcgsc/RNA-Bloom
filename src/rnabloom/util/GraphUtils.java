/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import rnabloom.RNABloom.ReadPair;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.Kmer;
import static rnabloom.util.SeqUtils.*;

/**
 *
 * @author gengar
 */
public final class GraphUtils {

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
        int range = end-start+1;
        
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
        return greedyExtendRightOnce(graph, source.getSuccessors(graph.getK(), graph.getMaxNumHash(), graph), lookahead, bf);
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
        return greedyExtendLeftOnce(graph, source.getPredecessors(graph.getK(), graph.getMaxNumHash(), graph), lookahead, bf);
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
    
    public static boolean represented(final ArrayList<Kmer> kmers,
                                    final BloomFilterDeBruijnGraph graph,
                                    final BloomFilter bf,
                                    final int lookahead,
                                    final int maxIndelSize,
                                    final int maxTipLength,
                                    final float percentIdentity) {
        int numKmers = kmers.size();
        int maxIndex = numKmers - 1;
        
        int k = graph.getK();
        int maxNumBubbleKmers = 3*k;
        
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
                
                if (assembledRange >= 3) {                    
                    if (startIndex > 0) {
                                                
                        if (lastRepresentedKmerFoundIndex < 0) {                            
                            if (startIndex >= maxTipLength || hasDepthLeft(kmers.get(0), graph, maxTipLength-startIndex)) {
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
                if (expectedLen >= maxTipLength || hasDepthRight(kmers.get(maxIndex), graph, maxTipLength-expectedLen)) {                    
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
        frontier.addAll(left.getSuccessors(k, numHash, graph));
        
        HashSet<String> kmersInFrontier = new HashSet<>();
        ArrayDeque<Kmer> newFrontier;
        for (int i=1; i<lowerBound; ++i) {
            kmersInFrontier.clear();
            newFrontier = new ArrayDeque<>();
            for (Kmer kmer : frontier) {
                if (bf.lookup(kmer.getHash())) {
                    for (Kmer s : kmer.getSuccessors(k, numHash, graph)) {
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
                if (bf.lookup(kmer.getHash())) {
                    if (kmer.equals(right)) {
                        return true;
                    }
                    newFrontier.add(kmer);
                    for (Kmer s : kmer.getSuccessors(k, numHash, graph)) {
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
            if (bf.lookup(kmer.getHash())) {
                if (kmer.equals(right)) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    public static ArrayList<Kmer> getSimilarCoveragePath(BloomFilterDeBruijnGraph graph, 
                                                        ArrayList<Kmer> leftKmers, 
                                                        ArrayList<Kmer> rightKmers, 
                                                        int bound, 
                                                        int lookahead, 
                                                        float maxCovGradient,
                                                        boolean rescue) {
        
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
            best.getSuccessors(k, numHash, graph, neighbors);

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
            best.getPredecessors(k, numHash, graph, neighbors);
            
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
    
    public static ArrayDeque<Kmer> findPath(BloomFilterDeBruijnGraph graph, Kmer left, Kmer right, int bound, int lookahead) {
        if (!graph.isLowComplexity(left) && !graph.isLowComplexity(right)) {
            int k = graph.getK();
            int numHash = graph.getMaxNumHash();

            ArrayDeque<Kmer> rightExtension = new ArrayDeque<>();

            ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
            for (int i=0; i<bound; ++i) {
                right.getPredecessors(k, numHash, graph, neighbors);

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
     * @return 
     */
    public static ArrayDeque<Kmer> getMaxCoveragePath(BloomFilterDeBruijnGraph graph, Kmer left, Kmer right, int bound, int lookahead) {
        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        HashSet<Kmer> leftPathKmers = new HashSet<>(bound);
        
        /* extend right */
        ArrayDeque<Kmer> leftPath = new ArrayDeque<>(bound);
        Kmer best;
        ArrayDeque<Kmer> neighbors;
        
        best = left;

        for (int depth=0; depth < bound; ++depth) {
            neighbors = best.getSuccessors(k, numHash, graph);
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
            neighbors = best.getPredecessors(k, numHash, graph);
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
                            
                            ArrayDeque<Kmer> path = getMaxCoveragePath(graph, kmers2.get(end+1), kmers.get(start-1), bubbleLength + maxIndelSize, lookahead);
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
        
        corrected = correctMismatches(kmers2, graph, covThreshold) || corrected;
        
        if (corrected) {
            return kmers2;
        }
        
        return null;
    }
    
    private static ArrayList<Kmer> correctErrorHelper(ArrayList<Kmer> kmers,
                                                    BloomFilterDeBruijnGraph graph, 
                                                    int lookahead,
                                                    int maxIndelSize,
                                                    float covThreshold,
                                                    float percentIdentity) {
        
        
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
                        ArrayDeque<Kmer> leftVars = kmers.get(i-1).getLeftVariants(k, numHash, graph);
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
                            m3 = getMinMedMaxKmerCoverage(testKmers);
                            
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
                            for (int j=i-numBadKmersSince; j<i; ++j) {
                                kmers2.add(kmers.get(j));
                            }                            
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
            ArrayDeque<Kmer> rightVars = kmers.get(i).getRightVariants(k, numHash, graph);
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
        
        corrected = correctMismatches(kmers2, graph, covThreshold) || corrected;
         
        if (corrected) {
            return kmers2;
        }
        
        // no changes
        return null;
    }
    
    public static boolean correctMismatches(ArrayList<Kmer> kmers,
                                            BloomFilterDeBruijnGraph graph,
                                            float covThreshold) {
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
                            float[] m = getMinMedMaxKmerCoverage(altKmers);
                            if (m[0] > 0 && m[1] > bestCov) {
                                bestCov = m[1];
                                bestAlt = altKmers;
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
                            float[] m = getMinMedMaxKmerCoverage(altKmers);
                            if (m[0] > 0 && m[1] > bestCov) {
                                bestCov = m[1];
                                bestAlt = altKmers;
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
                                    float percentIdentity) {
        
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
                                            percentIdentity);
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
            startIndex = numRightKmers - 1 - numFalsePositivesAllowed;
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
    
    public static ArrayDeque<ArrayList<Kmer>> breakWithReadPairedKmers(ArrayList<Kmer> kmers, BloomFilterDeBruijnGraph graph) {
        /**@TODO*/
        
        ArrayDeque<ArrayList<Kmer>> segments = new ArrayDeque<>();
        
        int d = graph.getReadKmerDistance();
        int lastIndex = kmers.size() - 1 - d;
        
        int start = -1;
        int end = -1;
        
        for (int i=0; i<=lastIndex; ++i) {
            if (graph.lookupReadKmerPair(kmers.get(i), kmers.get(i+d))) {
                if (start < 0) {
                    start = i;
                }
                
                end = i+d;
            }
            else if (start >= 0 && i >= end) {
                ArrayList<Kmer> sublist = new ArrayList<>(end-start+1);
                for (int j=start; j<=end; ++j) {
                    sublist.add(kmers.get(j));
                }
                segments.add(sublist);
                
                start = -1;
                end = -1;
            }
        }
        
        if (start >= 0) {
            ArrayList<Kmer> sublist = new ArrayList<>(end-start+1);
            for (int j=start; j<=end; ++j) {
                sublist.add(kmers.get(j));
            }
            segments.add(sublist);
        }
        
        return segments;
    }
    
    public static ArrayDeque<ArrayList<Kmer>> breakWithPairedKmers(ArrayList<Kmer> kmers, BloomFilterDeBruijnGraph graph) {
        /**@TODO Adjust how much paired kmers should interlock*/
        /**@TODO Adjust how many consecutive paired kmers are required*/
        
        ArrayDeque<ArrayList<Kmer>> segments = new ArrayDeque<>();
        
        int d = graph.getPairedKmerDistance();
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
            else if (start >= 0 && i >= end) {
                ArrayList<Kmer> sublist = new ArrayList<>(end-start+1);
                for (int j=start; j<=end; ++j) {
                    sublist.add(kmers.get(j));
                }
                segments.add(sublist);
                
                start = -1;
                end = -1;
            }
        }
        
        if (start >= 0) {
            ArrayList<Kmer> sublist = new ArrayList<>(end-start+1);
            for (int j=start; j<=end; ++j) {
                sublist.add(kmers.get(j));
            }
            segments.add(sublist);
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
        
        return leftWing + graph.assemble(pathKmers) + rightWing;
    }
        
    public static ArrayList<Kmer> overlap(ArrayList<Kmer> leftKmers, ArrayList<Kmer> rightKmers, BloomFilterDeBruijnGraph graph, int minOverlap) {
//        int k = graph.getK();
        
        String left = graph.assemble(leftKmers);
        String right = graph.assemble(rightKmers);
        String overlapped = overlapMaximally(left, right, minOverlap);
        
        if (overlapped != null) {
            int overlappedSeqLength = overlapped.length();
            int k = graph.getK();
                        
            int leftLen = left.length();
            int rightLen = right.length();
            
            if (overlappedSeqLength <= leftLen + rightLen - k) {
                // The overlap is larger than or equal to k
                
                boolean hasComplexKmer = false;
//                int end = leftLen - k + 1;
//                int start = overlappedKmers.size() - (rightLen - k + 1);

                int end = rightLen - (overlappedSeqLength - leftLen) - k + 1;

                for (int i=0; i<end; ++i) {
                    if (!graph.isLowComplexity(rightKmers.get(i))) {
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
            else {
                // The overlap is smaller than k
                                
                boolean hasComplexKmer = false;
                int start = leftLen - k + 1;
                int end = overlappedSeqLength - (rightLen - k + 1);

                ArrayList<Kmer> spanningKmers = graph.getKmers(overlapped, start, end);
                for (Kmer kmer : spanningKmers) {
                    if (kmer.count <= 0) {
                        // the overlap is not a valid path in the graph
                        return null;
                    }
                    
                    if (!hasComplexKmer && !graph.isLowComplexity(kmer)) {
                        // Require at least one complex kmers in the overlap
                        hasComplexKmer = true;
                    }
                }

                if (!hasComplexKmer) {
                    return null;
                }
                
                ArrayList<Kmer> overlappedKmers = new ArrayList<>(overlappedSeqLength - k + 1); //graph.getKmers(overlapped);
                overlappedKmers.addAll(leftKmers);
                overlappedKmers.addAll(spanningKmers);
                overlappedKmers.addAll(rightKmers);
                
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
                                                    boolean rescueSearch) {

        // 1. Attempt simple overlap
        ArrayList<Kmer> fragmentKmers = overlap(leftKmers, rightKmers, graph, minOverlap);
        
        if (fragmentKmers == null) {

            fragmentKmers = getSimilarCoveragePath(graph, leftKmers, rightKmers, bound, lookahead, maxCovGradient, rescueSearch);
//            ArrayDeque<Kmer> connectedPath = getMaxCoveragePath(graph, leftLastKmer, rightFirstKmer, bound, lookahead);

//            if (connectedPath == null && exhaustiveSearch) {
//                connectedPath = findPath(graph, leftLastKmer, rightFirstKmer, bound, lookahead);
//            }

//            if (connectedPath != null) {
//                fragmentKmers = new ArrayList<>(leftKmers.size() + connectedPath.size() + rightKmers.size());
//                fragmentKmers.addAll(leftKmers);
//                fragmentKmers.addAll(connectedPath);
//                fragmentKmers.addAll(rightKmers);
//            }
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
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        LinkedList<Kmer> results = new LinkedList<>();
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : source.getSuccessors(k, numHash, graph)) {
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
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        LinkedList<Kmer> results = new LinkedList<>();
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : source.getPredecessors(k, numHash, graph)) {
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
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        LinkedList<Kmer> results = new LinkedList<>();
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : source.getSuccessors(k, numHash, graph)) {
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
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        LinkedList<Kmer> results = new LinkedList<>();
        LinkedList<Float> values = new LinkedList<>();
        
        ListIterator<Kmer> resultsItr;
        ListIterator<Float> valuesItr;
        
        for (Kmer n : source.getPredecessors(k, numHash, graph)) {
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
        
        int distance = graph.getPairedKmerDistance();
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
        
        final int distance = graph.getPairedKmerDistance();
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
        
        final int distance = graph.getPairedKmerDistance();
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
                                            HashSet<Kmer> usedKmers) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        int distance = graph.getPairedKmerDistance();
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
                                            percentIdentity);
                        
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
        
        int distance = graph.getPairedKmerDistance();
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
                                            HashSet<Kmer> usedKmers) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        int distance = graph.getPairedKmerDistance();
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
                                            percentIdentity);
 
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
                                            float maxCovGradient) {
        
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        final int distance = graph.getPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int numKmers = kmers.size();
                            
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = cursor.getSuccessors(k, numHash, graph);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        int partnerIndex = numKmers-distance;
        final int minFragmentCoverageThreshold = (int) Math.floor(1.0/maxCovGradient);
        
        while (!neighbors.isEmpty()) {
            ArrayDeque<Kmer> simpleExtension = extendRight(cursor,
                                            graph, 
                                            maxTipLength, 
                                            usedKmers, 
                                            maxIndelSize, 
                                            percentIdentity);
                        
            if (!simpleExtension.isEmpty()) {
//                Iterator<Kmer> itr = simpleExtension.descendingIterator();
//                while(itr.hasNext() && assembledKmersBloomFilter.lookup(itr.next().getHash())) {
//                    itr.remove();
//                }
//                
//                if (!simpleExtension.isEmpty()) {
                    kmers.addAll(simpleExtension);
////                    usedKmers.addAll(simpleExtension);
                    
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
                    
                    float medEdgeCoverage = getMedianKmerCoverage(kmers, numKmers-1-lookahead, numKmers-1);
                    float minFragmentCoverage = getMinimumKmerCoverage(kmers, numKmers-distance, numKmers-1);
                    if (minFragmentCoverage >= minFragmentCoverageThreshold) {
                        ArrayDeque<Kmer> neighborsBackUp = new ArrayDeque<>(neighbors);
                        
                        // only remove neighbors based on coverage when fragment coverage is not too low
                        float minCovThreshold = minFragmentCoverage * maxCovGradient;            
                        itr = neighbors.iterator();
                        while (itr.hasNext()) {
                            Kmer kmer = itr.next();

                            if (kmer.count <= minCovThreshold) {
                                itr.remove();
                            }
                        }
                        
                        if (neighbors.size() > 1 && medEdgeCoverage > minFragmentCoverage) {
                            minCovThreshold = medEdgeCoverage * maxCovGradient;
                            
                            itr = neighbors.iterator();
                            while (itr.hasNext()) {
                                Kmer kmer = itr.next();

                                if (kmer.count <= minCovThreshold) {
                                    itr.remove();
                                }
                            }
                        }

                        if (neighbors.size() == 1) {
                            cursor = neighbors.pop();
                        }
                        else {
                            if (neighbors.isEmpty()) {
                                neighbors = neighborsBackUp;
                            }
                            
                            Kmer best = null;
                            float bestCov = 0;
                            float secondCov = 0;
                            
                            for (Kmer kmer : neighbors) {
                                
                                ArrayDeque<Kmer> e = extendRightWithPairedKmersOnly(kmer, kmers, graph, lookahead);
                                e.addFirst(kmer);
                                
                                float c = getMedianKmerCoverage(e);
                                
                                if (best == null) {
                                    best = kmer;
                                    bestCov = c;
                                } 
                                else {
                                    if (bestCov < c) {
                                        secondCov = bestCov;
                                        
                                        best = kmer;
                                        bestCov = c;
                                    }
                                    else if (secondCov < c) {
                                        secondCov = c;
                                    }
                                }
                            }

                            if (bestCov * maxCovGradient >= secondCov) {
                                cursor = best;
                            }
                            else {
                                return false; // ambiguous branches supported by paired kmers
                            }
                        }
                    }
                    else {
                        return false; // ambiguous branches supported by paired kmers
                    }
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
                                            float maxCovGradient) {
        final int k = graph.getK();
        final int numHash = graph.getMaxNumHash();
        
        final int distance = graph.getPairedKmerDistance();
//        int distanceInversePI = Math.max((int) (distance * (1-percentIdentity)), graph.getK());
        int numKmers = kmers.size();
        
        // Note that `kmers` are in reverse order already
        Kmer cursor = kmers.get(numKmers-1);
        ArrayDeque<Kmer> neighbors = cursor.getPredecessors(k, numHash, graph);
        if (neighbors.isEmpty()) {
            return false;
        }
        
        int partnerIndex = numKmers-distance;
        final int minFragmentCoverageThreshold = (int) Math.floor(1.0/maxCovGradient);
        
        while (!neighbors.isEmpty()) {
            ArrayDeque<Kmer> simpleExtension = extendLeft(cursor,
                                            graph, 
                                            maxTipLength, 
                                            usedKmers, 
                                            maxIndelSize, 
                                            percentIdentity);
            
            if (!simpleExtension.isEmpty()) {
//                Iterator<Kmer> itr = simpleExtension.descendingIterator();
//                while(itr.hasNext() && assembledKmersBloomFilter.lookup(itr.next().getHash())) {
//                    itr.remove();
//                }
//                
//                if (!simpleExtension.isEmpty()) {
                    kmers.addAll(simpleExtension);
//                    usedKmers.addAll(simpleExtension);
                    
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
                    
                    float medEdgeCoverage = getMedianKmerCoverage(kmers, numKmers-1-lookahead, numKmers-1);
                    float minFragmentCoverage = getMinimumKmerCoverage(kmers, numKmers-distance, numKmers-1);
                    if (minFragmentCoverage >= minFragmentCoverageThreshold) {
                        ArrayDeque<Kmer> neighborsBackUp = new ArrayDeque<>(neighbors);
                        
                        // only remove neighbors based on coverage when fragment coverage is not too low
                        float minCovThreshold = minFragmentCoverage * maxCovGradient;            
                        itr = neighbors.iterator();
                        while (itr.hasNext()) {
                            Kmer kmer = itr.next();

                            if (kmer.count <= minCovThreshold) {
                                itr.remove();
                            }
                        }
                        
                        if (neighbors.size() > 1 && medEdgeCoverage > minFragmentCoverage) {
                            minCovThreshold = medEdgeCoverage * maxCovGradient;
                            
                            itr = neighbors.iterator();
                            while (itr.hasNext()) {
                                Kmer kmer = itr.next();

                                if (kmer.count <= minCovThreshold) {
                                    itr.remove();
                                }
                            }
                        }
                        
                        if (neighbors.size() == 1) {
                            cursor = neighbors.pop();
                        }
                        else {
                            if (neighbors.isEmpty()) {
                                neighbors = neighborsBackUp;
                            }
                            
                            Kmer best = null;
                            float bestCov = 0;
                            float secondCov = 0;
                            
                            for (Kmer kmer : neighbors) {
                                
                                ArrayDeque<Kmer> e = extendLeftWithPairedKmersOnly(kmer, kmers, graph, lookahead);
                                e.addFirst(kmer);
                                
                                float c = getMedianKmerCoverage(e);
                                
                                if (best == null) {
                                    best = kmer;
                                    bestCov = c;
                                }
                                else {
                                    if (bestCov < c) {
                                        secondCov = bestCov;
                                        
                                        best = kmer;
                                        bestCov = c;
                                    }
                                    else if (secondCov < c) {
                                        secondCov = c;
                                    }
                                }
                            }

                            if (bestCov * maxCovGradient >= secondCov) {
                                cursor = best;
                            }
                            else {
                                return false; // ambiguous branches supported by paired kmers
                            }
                        }
                    }
                    else {
                        return false; // ambiguous branches supported by paired kmers
                    }
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
        
        final int distance = graph.getPairedKmerDistance();
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
        final int distance = graph.getPairedKmerDistance();
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
        final int distance = graph.getPairedKmerDistance();
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
        final int distance = graph.getPairedKmerDistance();
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
        final int distance = graph.getPairedKmerDistance();
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
        final int distance = graph.getPairedKmerDistance();
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
                                            int minNumPairs) {
        
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
                                    usedKmers);
        
        Collections.reverse(kmers);
        
        // extend with paired kmers RIGHT
        
        extendRightWithPairedKmersBFS(kmers, 
                                    graph,
                                    maxTipLength,
                                    maxIndelSize,
                                    percentIdentity,
                                    minNumPairs,
                                    usedKmers);
        
    }
    
    public static void extendWithPairedKmers(ArrayList<Kmer> kmers, 
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
            extendable = extendLeftWithPairedKmersBFS(kmers, 
                                        graph, 
                                        lookahead, 
                                        maxTipLength,
                                        maxIndelSize,
                                        percentIdentity,
                                        minNumPairs,
                                        assembledKmersBloomFilter,
                                        usedKmers,
                                        maxCovGradient);

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
                                        maxCovGradient);
            
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
    
    public static ArrayDeque<Kmer> naiveExtendRight(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<Kmer> terminators) {        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        HashSet<Kmer> usedKmers = new HashSet<>();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        kmer.getSuccessors(k, numHash, graph, neighbors);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : best.getLeftVariants(k, numHash, graph)) {
//                if (hasDepthLeft(s, graph, maxTipLength)) {
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
//                    if (hasDepthRight(n, graph, maxTipLength)) {
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
            
            best.getSuccessors(k, numHash, graph, neighbors);
        }
        
        return result;
    }
    
    public static ArrayList<Kmer> naiveExtend(ArrayList<Kmer> kmers, BloomFilterDeBruijnGraph graph, int maxTipLength) {
        HashSet<Kmer> usedKmers = new HashSet<>(kmers);
        
        ArrayDeque<Kmer> leftExtension = naiveExtendLeft(kmers.get(0), graph, maxTipLength, usedKmers);
        usedKmers.addAll(leftExtension);
        
        ArrayDeque<Kmer> rightExtension = naiveExtendRight(kmers.get(kmers.size()-1), graph, maxTipLength, usedKmers);
        
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
    
    public static ArrayDeque<Kmer> naiveExtendLeft(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<Kmer> terminators) {        
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        HashSet<Kmer> usedKmers = new HashSet<>();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        kmer.getPredecessors(k, numHash, graph, neighbors);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : best.getRightVariants(k, numHash, graph)) {
//                if (hasDepthRight(s, graph, maxTipLength)) {
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
//                    if (hasDepthLeft(n, graph, maxTipLength)) {
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
            
            best.getPredecessors(k, numHash, graph, neighbors);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> extendRight(Kmer source,
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLength, 
                                            HashSet<Kmer> usedKmers, 
                                            int maxIndelSize, 
                                            float percentIdentity) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        source.getSuccessors(k, numHash, graph, neighbors);
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
                
                ArrayDeque<Kmer> b = naiveExtendRight(best, graph, maxTipLength, usedKmers);
                
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
                        ArrayDeque<Kmer> b = naiveExtendRight(n, graph, maxTipLength, usedKmers);
                        
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
                    best.getSuccessors(k, numHash, graph, neighbors);
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
            
            best.getSuccessors(k, numHash, graph, neighbors);
        }
        
        return result;
    }
    
    public static ArrayDeque<Kmer> extendLeft(Kmer source,
                                            BloomFilterDeBruijnGraph graph, 
                                            int maxTipLength, 
                                            HashSet<Kmer> usedKmers, 
                                            int maxIndelSize, 
                                            float percentIdentity) {
        int k = graph.getK();
        int numHash = graph.getMaxNumHash();
        
        ArrayDeque<Kmer> result = new ArrayDeque<>();
        
        ArrayDeque<Kmer> neighbors = new ArrayDeque<>(4);
        source.getPredecessors(k, numHash, graph, neighbors);
        Kmer best;
        
        while (!neighbors.isEmpty()) {
            
            if (neighbors.size() == 1) {
                best = neighbors.pop();
                
                if (usedKmers.contains(best)) {
                    return result;
                }
                
                result.add(best);
                
                ArrayDeque<Kmer> b = naiveExtendLeft(best, graph, maxTipLength, usedKmers);
                
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
                        ArrayDeque<Kmer> b = naiveExtendLeft(n, graph, maxTipLength, usedKmers);
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
                    best.getPredecessors(k, numHash, graph, neighbors);
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
            
            best.getPredecessors(k, numHash, graph, neighbors);
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
    
    public static boolean isFusion(ArrayList<Kmer> seqKmers, BloomFilterDeBruijnGraph graph, BloomFilter assembledKmers, int lookahead) {
        int k = graph.getK();   
        int numKmers = seqKmers.size();
        
        if (assembledKmers.lookup(seqKmers.get(0).getHash()) &&
                assembledKmers.lookup(seqKmers.get(numKmers-1).getHash())) {
            
            int i = 1;
            for (; i < numKmers-1; ++i) {
                if (!assembledKmers.lookup(seqKmers.get(i).getHash())) {
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
                    break;
                }
            }
            
            ++j;
            
            if (j - i <= k) {
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
    
    public static boolean isTemplateSwitch(ArrayList<Kmer> seqKmers, BloomFilterDeBruijnGraph graph, BloomFilter assembledKmers, int lookahead) {
        int k = graph.getK();   
        int numKmers = seqKmers.size();
        
        float leftEdgeCov = getMinimumKmerCoverage(seqKmers, 0, Math.min(lookahead, numKmers));
        float rightEdgeCov = getMinimumKmerCoverage(seqKmers, Math.max(0, numKmers-lookahead), numKmers);
        
        if (assembledKmers.lookup(seqKmers.get(0).getHash()) &&
                (!assembledKmers.lookup(seqKmers.get(numKmers-1).getHash()) || leftEdgeCov > rightEdgeCov)) {
            int i = 1;
            for (; i < numKmers-1; ++i) {
                if (!assembledKmers.lookup(seqKmers.get(i).getHash())) {
                    break;
                }
            }
            
            if (i == numKmers -1) {
                return false;
            }
            
            --i;
            
            if (i+k < numKmers) {
                String tipRC = reverseComplement(graph.assemble(seqKmers, i+k, numKmers));
                
                ArrayDeque<Kmer> leftExtension = greedyExtendLeft(graph, seqKmers.get(0), lookahead, 1000, assembledKmers);
                ArrayDeque<Kmer> rightExtension = greedyExtendRight(graph, seqKmers.get(i), lookahead, 1000, assembledKmers);
                
                leftExtension.addAll(seqKmers.subList(0, i+1));
                leftExtension.addAll(rightExtension);
                
                String backbone = graph.assemble(leftExtension);
                if (backbone.contains(tipRC)) {
                    return true;
                }
            }
        }
        else if (assembledKmers.lookup(seqKmers.get(numKmers-1).getHash()) &&
                (!assembledKmers.lookup(seqKmers.get(0).getHash()) || leftEdgeCov < rightEdgeCov)) {
            int j = numKmers-2;
            for (; j > 0; --j) {
                if (!assembledKmers.lookup(seqKmers.get(j).getHash())) {
                    break;
                }
            }
            
            ++j;
            
            if (j-k > 0) {
                String tipRC = reverseComplement(graph.assemble(seqKmers, 0, j-k));
                
                ArrayDeque<Kmer> leftExtension = greedyExtendLeft(graph, seqKmers.get(j), lookahead, 1000, assembledKmers);
                ArrayDeque<Kmer> rightExtension = greedyExtendRight(graph, seqKmers.get(numKmers-1), lookahead, 1000, assembledKmers);
                
                leftExtension.addAll(seqKmers.subList(j, numKmers));
                leftExtension.addAll(rightExtension);
                
                String backbone = graph.assemble(leftExtension);
                if (backbone.contains(tipRC)) {
                    return true;
                }
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
