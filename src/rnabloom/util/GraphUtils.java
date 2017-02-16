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
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import rnabloom.RNABloom.ReadPair;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.hash.NTHashIterator;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
import static rnabloom.util.SeqUtils.getFirstKmer;
import static rnabloom.util.SeqUtils.getLastKmer;
import static rnabloom.util.SeqUtils.getNumGC;
import static rnabloom.util.SeqUtils.kmerizeToCollection;
import static rnabloom.util.SeqUtils.overlapMaximally;

/**
 *
 * @author gengar
 */
public final class GraphUtils {
    
    public static float getMedianKmerCoverage(final ArrayList<Kmer> kmers) {
        int numKmers = kmers.size();
        int halfNumKmers = numKmers/2;
        
        ArrayList<Float> counts = new ArrayList<>(numKmers);
        for (Kmer kmer : kmers) {
            counts.add(kmer.count);
        }
        
        Collections.sort(counts);
        
        if (numKmers % 2 == 0) {
            return (counts.get(halfNumKmers) + counts.get(halfNumKmers -1))/2.0f;
        }
        
        return counts.get(halfNumKmers);
    }
    
    public static float getMaxMedianCoverageRight(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead) {
        Iterator<Kmer> itr = graph.getSuccessors(source).iterator();
        
        if (!itr.hasNext()) {
            return source.count;
        }
        else {
            Kmer cursor = itr.next();
            ArrayList<Kmer> path = new ArrayList<>(lookahead); 
            path.add(source);
            path.add(cursor);

            ArrayList<Iterator<Kmer>> frontier = new ArrayList<>(lookahead);
            frontier.add(itr);

            float bestCov = 0;

            while (!frontier.isEmpty()) {
                if (path.size() < lookahead) {
                    itr = graph.getSuccessors(cursor).iterator();
                    if (itr.hasNext()) {
                        cursor = itr.next();
                        path.add(cursor);
                        frontier.add(itr);
                        continue;
                    }
                }

                float pathCov = getMedianKmerCoverage(path);
                if (bestCov < pathCov) {
                    bestCov = pathCov;
                }

                int i = path.size()-2;
                while (i >= 0) {
                    itr = frontier.get(i);
                    path.remove(i+1);
                    if (!itr.hasNext()) {
                        frontier.remove(i);
                        --i;
                    }
                    else {
                        cursor = itr.next();
                        path.add(cursor);
                        break;
                    }
                }
            }
            
            return bestCov;
        }
    }

    public static float getMaxMedianCoverageLeft(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead) {
        Iterator<Kmer> itr = graph.getPredecessors(source).iterator();
        
        if (!itr.hasNext()) {
            return source.count;
        }
        else {
            Kmer cursor = itr.next();
            ArrayList<Kmer> path = new ArrayList<>(lookahead);
            path.add(source);
            path.add(cursor);

            ArrayList<Iterator<Kmer>> frontier = new ArrayList<>(lookahead);
            frontier.add(itr);

            float bestCov = 0;

            while (!frontier.isEmpty()) {
                if (path.size() < lookahead) {
                    itr = graph.getPredecessors(cursor).iterator();
                    if (itr.hasNext()) {
                        cursor = itr.next();
                        path.add(cursor);
                        frontier.add(itr);
                        continue;
                    }
                }

                float pathCov = getMedianKmerCoverage(path);
                if (bestCov < pathCov) {
                    bestCov = pathCov;
                }

                int i = path.size()-2;
                while (i >= 0) {
                    itr = frontier.get(i);
                    path.remove(i+1);
                    if (!itr.hasNext()) {
                        frontier.remove(i);
                        --i;
                    }
                    else {
                        cursor = itr.next();
                        path.add(cursor);
                        break;
                    }
                }
            }

            return bestCov;
        }
    }
    
    private static Kmer greedyExtendRightOnce(BloomFilterDeBruijnGraph graph, ArrayDeque<Kmer> candidates, int lookahead) {
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
    
    public static Kmer greedyExtendRightOnce(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead) {
        ArrayDeque<Kmer> candidates = graph.getSuccessors(source);
        
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
    
    public static Kmer greedyExtendLeftOnce(BloomFilterDeBruijnGraph graph, ArrayDeque<Kmer> candidates, int lookahead) {
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
    
    public static Kmer greedyExtendLeftOnce(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead) {
        ArrayDeque<Kmer> candidates = graph.getPredecessors(source);
        
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
    
    public static boolean hasValidPath(BloomFilterDeBruijnGraph graph, Kmer left, Kmer right, BloomFilter bf, int lowerBound, int upperBound) {
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
                        if (!kmersInFrontier.contains(s.seq)) { 
                            newFrontier.add(s);
                            kmersInFrontier.add(s.seq);
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
                        if (!kmersInFrontier.contains(s.seq)) { 
                            newFrontier.add(s);
                            kmersInFrontier.add(s.seq);
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
    public static ArrayList<Kmer> getMaxCoveragePath(BloomFilterDeBruijnGraph graph, Kmer left, Kmer right, int bound, int lookahead) {
        
        HashSet<String> leftPathKmers = new HashSet<>(bound);
        
        /* extend right */
        ArrayList<Kmer> leftPath = new ArrayList<>(bound);
        Kmer best = left;
        ArrayDeque<Kmer> neighbors;
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
                    best = greedyExtendRightOnce(graph, best, lookahead);
                }
                
                if (best.equals(right)) {
                    return leftPath;
                }
                else {
                    leftPath.add(best);
                }
            }
        }
        
        for (Kmer kmer : leftPath) {
            leftPathKmers.add(kmer.seq);
        }
        
        /* not connected, search from right */
        ArrayList<Kmer> rightPath = new ArrayList<>(bound);
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
                    best = greedyExtendLeftOnce(graph, best, lookahead);
                }
                
                if (best.equals(left)) {
                    Collections.reverse(rightPath);
                    return rightPath;
                }
                else if (leftPathKmers.contains(best.seq)) {
                    /* right path intersects the left path */
                    String convergingKmer = best.seq;
                    ArrayList<Kmer> path = new ArrayList<>(bound);
                    for (Kmer kmer : leftPath) {
                        path.add(kmer);
                        if (convergingKmer.equals(kmer.seq)) {
                            break;
                        }
                    }
                    
                    if (path.size() + rightPath.size() <= bound) {
                        Collections.reverse(rightPath);
                        path.addAll(rightPath);
                        return path;
                    }
                }
                else {
                    rightPath.add(best);
                }
            }
        }
        
        return null;
    }
    
    public static float getMedianModifyInput(float[] arr) {
        Arrays.sort(arr);
        int halfLen = arr.length/2;
        if (halfLen % 2 == 0) {
            return (arr[halfLen-1] + arr[halfLen])/2.0f;
        }
        
        return arr[halfLen];
    }    
    
    public static float getMedian(float[] arr) {
        float[] a = Arrays.copyOf(arr, arr.length);
        Arrays.sort(a);
        int halfLen = a.length/2;
        if (halfLen % 2 == 0) {
            return (a[halfLen-1] + a[halfLen])/2.0f;
        }
        
        return a[halfLen];
    }
        
    public static float getMedian(ArrayDeque<Float> arr) {
        ArrayList<Float> a = new ArrayList<>(arr);
        Collections.sort(a);
        int halfLen = a.size()/2;
        if (halfLen % 2 == 0) {
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
    
    public static ArrayList<Kmer> greedyExtendLeft(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead, int bound) {
        ArrayList<Kmer> extension = new ArrayList<>(bound);
        Kmer nextKmer = source;
        extension.add(nextKmer);
        
        for (int i=0; i<bound; ++i) {
            nextKmer = greedyExtendLeftOnce(graph, nextKmer, lookahead);
            extension.add(nextKmer);
            
            if (nextKmer == null) {
                break;
            }
        }
        
        Collections.reverse(extension);
        
        return extension;
    }
    
    public static ArrayList<Kmer> greedyExtendRight(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead, int bound) {
        ArrayList<Kmer> extension = new ArrayList<>(bound);
        Kmer nextKmer = source;
        extension.add(nextKmer);
        
        for (int i=0; i<bound; ++i) {
            nextKmer = greedyExtendRightOnce(graph, nextKmer, lookahead);
            extension.add(nextKmer);
            
            if (nextKmer == null) {
                break;
            }
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
    
    public static ReadPair correctErrors2(ArrayList<Kmer> leftKmers, 
                                            ArrayList<Kmer> rightKmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead,
                                            int maxIndelSize, 
                                            float maxCovGradient, 
                                            float covFPR,
                                            int errorCorrectionIterations) {
        
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

            // find cov threshold in left kmers
            boolean leftThresholdFound = false;
            int startIndex = numLeftKmers - 1 - numFalsePositivesAllowed;
            float leftCovThreshold = covs[startIndex];
            float c = -1;
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

            // find cov threshold in right kmers
            boolean rightThresholdFound = false;
            startIndex = numRightKmers - 1 - numFalsePositivesAllowed;
            float rightCovThreshold = covs[startIndex];
            c = -1;
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

                // correct left read
                leftCorrected = false;
                ArrayList<Kmer> leftKmers2 = new ArrayList<>(numLeftKmers);
                int numBadKmersSince = 0;
                Kmer kmer;
                for (int i=0; i<numLeftKmers; ++i) {
                    kmer = leftKmers.get(i);
                    if (kmer.count >= covThreshold) {                    
                        if (numBadKmersSince > 0) {
                            if (!leftKmers2.isEmpty()) {
                                ArrayList<Kmer> path = getMaxCoveragePath(graph, leftKmers2.get(leftKmers2.size()-1), kmer, numBadKmersSince + maxIndelSize, lookahead);
                                if (path == null) {
                                    for (int j=i-numBadKmersSince; j<i; ++j) {
                                        leftKmers2.add(leftKmers.get(j));
                                    }
                                }
                                else {
                                    leftKmers2.addAll(path);
                                    leftCorrected = true;
                                }
                            }

                            numBadKmersSince = 0;
                        }

                        leftKmers2.add(kmer);
                    }
                    else {
                        ++numBadKmersSince;
                    }
                }

                if (leftCorrected || leftKmers2.size() != numLeftKmers) {
                    leftKmers = leftKmers2;
                    leftCorrected = true;
                }
                
                // correct right read
                rightCorrected = false;
                ArrayList<Kmer> rightKmers2 = new ArrayList<>(numRightKmers);
                numBadKmersSince = 0;
                for (int i=0; i<numRightKmers; ++i) {
                    kmer = rightKmers.get(i);
                    if (kmer.count >= covThreshold) {                    
                        if (numBadKmersSince > 0) {
                            if (!rightKmers2.isEmpty()) {
                                ArrayList<Kmer> path = getMaxCoveragePath(graph, rightKmers2.get(rightKmers2.size()-1), kmer, numBadKmersSince + maxIndelSize, lookahead);
                                if (path == null) {
                                    for (int j=i-numBadKmersSince; j<i; ++j) {
                                        rightKmers2.add(rightKmers.get(j));
                                    }
                                }
                                else {
                                    rightKmers2.addAll(path);
                                    rightCorrected = true;
                                }
                            }

                            numBadKmersSince = 0;
                        }

                        rightKmers2.add(kmer);
                    }
                    else {
                        ++numBadKmersSince;
                    }
                }

                if (rightCorrected || rightKmers2.size() != numRightKmers) {
                    rightKmers = rightKmers2;
                    rightCorrected = true;
                }
                
                if (!leftCorrected && !rightCorrected) {
                    break;
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

    public static String assemble(ArrayList<Kmer> kmers, int k) {
        StringBuilder sb = new StringBuilder(kmers.size() + k - 1);

        int lastIndex = k - 1;        
        sb.append(kmers.get(0).seq.substring(0, lastIndex));        
        for (Kmer kmer : kmers) {
            sb.append(kmer.seq.charAt(lastIndex));
        }

        // surprisingly slower!
//        int numKmers = kmers.size();
//        int length = numKmers + k - 1;
//        StringBuilder sb = new StringBuilder(length);
//        int numNonOverlappingKmers = length/k;
//
//        for (int i=0; i<numNonOverlappingKmers; ++i) {
//            sb.append(kmers.get(i*k).seq);
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
    
    public static String assembleFirstBase(ArrayList<Kmer> kmers) {
        StringBuilder sb = new StringBuilder(kmers.size());
        for (Kmer kmer : kmers) {
            sb.append(kmer.seq.charAt(0));
        }
        
        return sb.toString();
    }

    public static String assembleLastBase(ArrayList<Kmer> kmers) {
        int lastIndex = kmers.get(0).seq.length() - 1;
        
        StringBuilder sb = new StringBuilder(kmers.size());
        for (Kmer kmer : kmers) {
            sb.append(kmer.seq.charAt(lastIndex));
        }
        
        return sb.toString();
    }
    
    public static ArrayList<Kmer> greedyExtend(Kmer seed, BloomFilterDeBruijnGraph graph, int lookahead) {
        /**@TODO store smallest strand kmers for non-strand specific sequences */
        
        HashSet<String> pathKmerStr = new HashSet<>(1000);
        pathKmerStr.add(seed.seq);
        
        ArrayList<Kmer> rightPath = new ArrayList<>(1000);
        
        /* extend on right side */
        Kmer best = seed;
        while (true) {
            best = greedyExtendRightOnce(graph, best, lookahead);
            if (best != null) {
                String seq = best.seq;
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
                String seq = best.seq;
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
    
    public static float getGCContent(ArrayList<Kmer> path, int k) {
        
        int pathLen = path.size();
        
        int numGC = getNumGC(path.get(0).seq);
        
        char c;
        int lastCharIndex = k-1;
        for (int i=1; i<pathLen; ++i) {
            c = path.get(i).seq.charAt(lastCharIndex);
            if (c == 'G' || c == 'C') {
                ++numGC;
            }
        }
        
        return (float) numGC / (k + pathLen -1);
    }
    
    public static ArrayList<Kmer> getKmers(String seq, BloomFilterDeBruijnGraph graph, int maxIndelSize, int lookahead) {
        int k = graph.getK();
        int numKmers = seq.length() - k + 1;
        
        ArrayList<Kmer> result = new ArrayList<>(numKmers);
        ArrayList<Kmer> bestResult = null;
        
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
                nextKmer = new Kmer(seq.substring(i, i+k), c, hVals);
                
                if (numMissingKmers > 0 && !result.isEmpty()) {
                    ArrayList<Kmer> path = getMaxCoveragePath(graph, result.get(result.size()-1), nextKmer, maxIndelSize + numMissingKmers, lookahead);
                    
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
                        ArrayList<Kmer> path = getMaxCoveragePath(graph, current.get(current.size()-1), next.get(0), bound, lookahead);

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
        
        ArrayList<Kmer> pathKmers = getMaxCoveragePath(graph, graph.getKmer(getLastKmer(left, k)), graph.getKmer(getFirstKmer(right, k)), bound, lookahead);
        
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
        
        String left = assemble(leftKmers, k);
        String right = assemble(rightKmers, k);
        
        String overlapped = overlapMaximally(left, right, minOverlap);
        
        if (overlapped != null) {
            int overlappedLength = overlapped.length();
            if (overlappedLength >= left.length() && overlappedLength >= right.length()) {
                ArrayList<Kmer> overlappedKmers = graph.getKmers(overlapped);
                for (Kmer kmer : overlappedKmers) {
                    if (kmer.count <= 0) {
                        return null;
                    }
                }

                return overlappedKmers;
            }
        }
        
        return null;
    }
    
    public static String overlapThenConnect(String left, String right, BloomFilterDeBruijnGraph graph, int bound, int lookahead, int minOverlap) {
        
        // overlap before finding path
        String overlapped = overlapMaximally(left, right, minOverlap);
        if (overlapped != null && graph.isValidSeq(overlapped)) {
            return overlapped;
        }
        
        return connect(left, right, graph, bound, lookahead);
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
    
    private static int maxRightPartnerSearchDepth(ArrayList<Kmer> fragmentKmers, BloomFilterDeBruijnGraph graph, int pairedKmerDistance) {
        
        final int numKmers = fragmentKmers.size();
       
        for (int i=numKmers-graph.getK(); i>Math.max(0, numKmers-pairedKmerDistance); --i) {
            if (graph.lookupLeftKmer(fragmentKmers.get(i).hashVals)) {
                if (i>0 && graph.lookupLeftKmer(fragmentKmers.get(i-1).hashVals)) {
                    return pairedKmerDistance - (numKmers - i);
                }
            }
        }
        
        return 0;
    }
    
    private static int maxLeftPartnerSearchDepth(ArrayList<Kmer> fragmentKmers, BloomFilterDeBruijnGraph graph, int pairedKmerDistance) {
        
        final int numKmers = fragmentKmers.size();
       
        for (int i=numKmers-graph.getK(); i>Math.max(0, numKmers-pairedKmerDistance); --i) {
            if (graph.lookupRightKmer(fragmentKmers.get(i).hashVals)) {
                if (i>0 && graph.lookupLeftKmer(fragmentKmers.get(i-1).hashVals)) {
                    return pairedKmerDistance - (numKmers - i);
                }
            }
        }
        
        return 0;
    }
    
    public static void extendWithPairedKmers(ArrayList<Kmer> kmers, 
                                            BloomFilterDeBruijnGraph graph, 
                                            int lookahead, 
                                            int maxTipLength, 
                                            boolean greedy, 
                                            BloomFilter assembledKmersBloomFilter) {
        
        final int distance = graph.getPairedKmerDistance();
        
        int maxRightDepth = maxRightPartnerSearchDepth(kmers, graph, distance);
                
        // kmer pairs used in the extension of this transcript
        final HashSet<String> usedPairs = new HashSet<>();
        
        final HashSet<String> usedKmers = new HashSet<>();
        for (Kmer kmer : kmers) {
            usedKmers.add(kmer.seq);
        }
        
        // extend right
  
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

        while (!branchesStack.isEmpty() && !stop && maxRightDepth > 0) {
            neighbors = branchesStack.getLast();
            depth = extension.size();
            
            if (neighbors.isEmpty()) {
                branchesStack.removeLast();
//                    extension.removeLast().successors = null; // prune the cached successors
                extension.pollLast();
            }
            else if (depth == 0 && neighbors.size() == 1 && !usedKmers.contains(neighbors.peek().seq)) {
                // naive extension                
                n = neighbors.removeFirst();
                kmers.add(n);
                usedKmers.add(n.seq);
                branchesStack.add(getSuccessorsRanked(n, graph, lookahead));
//                n.successors = null; // prune the cached successors
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

                                    kmers.addAll(extension);
                                    kmers.add(n);
                                    
                                    if (!greedy && 
                                            assembledKmersBloomFilter.lookup(n.hashVals) &&
                                            assembledKmersBloomFilter.lookup(partner.hashVals) &&
                                            assembledKmersBloomFilter.lookup(partner2.hashVals) &&
                                            assembledKmersBloomFilter.lookup(n2.hashVals) ) {
                                        stop = true;
                                        break;
                                    }
                                    
                                    Iterator<Kmer> kmerItr = extension.iterator();
                                    while (kmerItr.hasNext()) {
                                        kmer = kmerItr.next();
        //                                kmer.successors = null; // prune the cached successors
                                        usedKmers.add(kmer.seq);
                                        kmerItr.remove();
                                    }
                                    usedKmers.add(n.seq);

                                    visitedKmers.clear();

                                    mergedSeq = partner.seq + n.seq;

                                    if (usedPairs.contains(mergedSeq)) {
                                        stop = true;
                                        break;
                                    }
                                    else {
                                        usedPairs.add(mergedSeq);
                                    }

                                    maxRightDepth = maxRightPartnerSearchDepth(kmers, graph, distance);
                                    branchesStack.add(getSuccessorsRanked(n, graph, lookahead));
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
                            
                            kmers.addAll(extension);
                            kmers.add(n);
                            
                            if (!greedy && 
                                    assembledKmersBloomFilter.lookup(n.hashVals) &&
                                    assembledKmersBloomFilter.lookup(partner.hashVals) &&
                                    assembledKmersBloomFilter.lookup(partner2.hashVals) &&
                                    assembledKmersBloomFilter.lookup(n2.hashVals)) {
                                stop = true;
                                break;
                            }
                            
                            Iterator<Kmer> kmerItr = extension.iterator();
                            while (kmerItr.hasNext()) {
                                kmer = kmerItr.next();
//                                kmer.successors = null; // prune the cached successors
                                usedKmers.add(kmer.seq);
                                kmerItr.remove();
                            }
                            usedKmers.add(n.seq);
                            
                            visitedKmers.clear();

                            if (usedPairs.contains(mergedSeq)) {
                                stop = true;
                                break;
                            }
                            else {
                                usedPairs.add(mergedSeq);
                            }
                            
                            maxRightDepth = maxRightPartnerSearchDepth(kmers, graph, distance);
                            branchesStack.add(getSuccessorsRanked(n, graph, lookahead));
                            
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
        
        // extend left
        
        Collections.reverse(kmers);
        
        stop = false;
        visitedKmers.clear();
        extension.clear();
        branchesStack.clear();
        neighbors = getPredecessorsRanked(kmers.get(kmers.size()-1), graph, lookahead);
        branchesStack.add(neighbors);
        int maxLeftDepth = maxLeftPartnerSearchDepth(kmers, graph, distance);
        
        while (!branchesStack.isEmpty() && !stop && maxLeftDepth > 0) {
            neighbors = branchesStack.getLast();
            depth = extension.size();
            
            if (neighbors.isEmpty()) {
                branchesStack.removeLast();
//                    extension.removeLast().predecessors = null; // prune the cached predecessors
                extension.pollLast();
            }
            else if (depth == 0 && neighbors.size() == 1 && !usedKmers.contains(neighbors.peek().seq)) {
                // naive extension
                n = neighbors.removeFirst();
                kmers.add(n);
                usedKmers.add(n.seq);
                branchesStack.add(getPredecessorsRanked(n, graph, lookahead));
//                n.predecessors = null; // prune the cached predecessors
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

                                    if (!greedy && 
                                            assembledKmersBloomFilter.lookup(n.hashVals) && 
                                            assembledKmersBloomFilter.lookup(partner.hashVals)) {
                                        stop = true;
                                        break;
                                    }

                                    branchesStack.clear();

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

                                    visitedKmers.clear();

                                    mergedSeq = n.seq + partner.seq;

                                    if (usedPairs.contains(mergedSeq)) {
                                        stop = true;
                                        break;
                                    }
                                    else {
                                        usedPairs.add(mergedSeq);
                                    }

                                    maxLeftDepth = maxLeftPartnerSearchDepth(kmers, graph, distance);
                                    branchesStack.add(getPredecessorsRanked(n, graph, lookahead));

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
                            
                            if (!greedy && 
                                    assembledKmersBloomFilter.lookup(n.hashVals) && 
                                    assembledKmersBloomFilter.lookup(partner.hashVals)) {
                                stop = true;
                                break;
                            }

                            branchesStack.clear();

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
                            
                            visitedKmers.clear();

                            if (usedPairs.contains(mergedSeq)) {
                                stop = true;
                            }
                            else {
                                usedPairs.add(mergedSeq);
                            }

                            maxLeftDepth = maxLeftPartnerSearchDepth(kmers, graph, distance);
                            branchesStack.add(getPredecessorsRanked(n, graph, lookahead));
                            
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
        
        Collections.reverse(kmers);
        
        //return assemble(kmers, k);
    }
    
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
    
    private static boolean hasDepthRight(String source, BloomFilterDeBruijnGraph graph, int depth) {
        ArrayDeque<ArrayDeque> frontier = new ArrayDeque<>();
        ArrayDeque<String> alts = graph.getSuccessors(source);
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
    
    private static boolean hasDepthRight(Kmer source, BloomFilterDeBruijnGraph graph, int depth) {
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
    
    private static boolean hasDepthLeft(String source, BloomFilterDeBruijnGraph graph, int depth) {
        ArrayDeque<ArrayDeque> frontier = new ArrayDeque<>();
        ArrayDeque<String> alts = graph.getPredecessors(source);
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
    
    private static boolean hasDepthLeft(Kmer source, BloomFilterDeBruijnGraph graph, int depth) {
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
    
    public static String naiveExtendWithSimpleBubblePopping(String fragment, BloomFilterDeBruijnGraph graph, int maxTipLength) {
        int k = graph.getK();
        int kMinus1 = k - 1;
        HashSet<String> fragmentKmers = new HashSet<>(2*(fragment.length()-k+1));
        kmerizeToCollection(fragment, k, fragmentKmers);
        
        ArrayDeque<String> neighbors = graph.getPredecessors(getFirstKmer(fragment, k));
        while (neighbors.size() > 0) {
            String bestExtension = null;
            String bestPrefix = null;
            float bestCov = 0;
            
            for (String p : neighbors) {
                String e = naiveExtendLeft(p, graph, maxTipLength, fragmentKmers);
                if (e.length() > maxTipLength) {
                    String ext = e + p;
                    
                    if (bestExtension == null) {
                        bestExtension = e + p.charAt(0);
                        bestPrefix = getFirstKmer(ext, kMinus1);
                        bestCov = getMedianKmerCoverage(graph.getKmers(ext));
                    }
                    else {
                        String myPrefix = getFirstKmer(ext, kMinus1);
                        if (bestPrefix.equals(myPrefix)) {
                            if (Math.abs(e.length() + 1 - bestExtension.length()) <= 1) {
                                float myCov = getMedianKmerCoverage(graph.getKmers(ext));
                                if (myCov > bestCov) {
                                    bestExtension = e + p.charAt(0);
                                    bestPrefix = myPrefix;
                                    bestCov = myCov;
                                }
                            }
                            else {
                                bestExtension = null;
                                bestPrefix = null;
                                bestCov = 0;
                                break;
                            }
                        }
                        else {
                            bestExtension = null;
                            bestPrefix = null;
                            bestCov = 0;
                            break;
                        }
                    }
                }
            }
            
            if (bestExtension == null) {
                break;
            }
            else {
                fragment = bestExtension + fragment;
            }
            
            /**@TODO extend to junction kmer*/
            
            neighbors = graph.getPredecessors(getFirstKmer(fragment, k));
        }
        
        neighbors = graph.getSuccessors(getLastKmer(fragment, k));
        while (neighbors.size() > 0) {
            String bestExtension = null;
            String bestPrefix = null;
            float bestCov = 0;
            
            for (String s : neighbors) {
                String e = naiveExtendRight(s, graph, maxTipLength, fragmentKmers);
                if (e.length() > maxTipLength) {
                    String ext = s + e;
                    
                    if (bestExtension == null) {
                        bestExtension = s.charAt(kMinus1) + e;
                        bestPrefix = getLastKmer(ext, kMinus1);
                        bestCov = getMedianKmerCoverage(graph.getKmers(ext));
                    }
                    else {
                        String myPrefix = getLastKmer(ext, kMinus1);
                        if (bestPrefix.equals(myPrefix)) {
                            if (Math.abs(e.length() + 1 - bestExtension.length()) <= 1) {
                                float myCov = getMedianKmerCoverage(graph.getKmers(ext));
                                if (myCov > bestCov) {
                                    bestExtension = s.charAt(kMinus1) + e;
                                    bestPrefix = myPrefix;
                                    bestCov = myCov;
                                }
                            }
                            else {
                                bestExtension = null;
                                bestPrefix = null;
                                bestCov = 0;
                                break;
                            }
                        }
                        else {
                            bestExtension = null;
                            bestPrefix = null;
                            bestCov = 0;
                            break;
                        }
                    }
                }
            }
            
            if (bestExtension == null) {
                break;
            }
            else {
                fragment = fragment + bestExtension;
            }
            
            /**@TODO extend to junction kmer*/
            
            neighbors = graph.getSuccessors(getLastKmer(fragment, k));
        }
        
        return fragment;
    }
    
    public static String naiveExtend(String fragment, BloomFilterDeBruijnGraph graph, int maxTipLength) {
        int k = graph.getK();
        HashSet<String> fragmentKmers = new HashSet<>(2*(fragment.length()-k+1));
        kmerizeToCollection(fragment, k, fragmentKmers);
        
        return naiveExtendLeft(getFirstKmer(fragment, k), graph, maxTipLength, fragmentKmers) + fragment + naiveExtendRight(getLastKmer(fragment, k), graph, maxTipLength, fragmentKmers);
    }
    
    private static String naiveExtendRight(String kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<String> terminators) {        
        StringBuilder sb = new StringBuilder(100);
        int lastBaseIndex = graph.getKMinus1();
        
        ArrayDeque<String> neighbors = graph.getSuccessors(kmer);
        String best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (String s : graph.getLeftVariants(best)) {
                if (hasDepthLeft(s, graph, maxTipLength)) {
                    return sb.toString();
                }
            }
            
            if (neighbors.size() == 1) {
                best = neighbors.peek();
            }
            else {
                best = null;
                for (String n : neighbors) {
                    if (hasDepthRight(n, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            return sb.toString();
                        }
                    }
                }
            }
            
            if (best == null || terminators.contains(best)) {
                break;
            }
            
            sb.append(best.charAt(lastBaseIndex));
            terminators.add(best);
        }
        
        return sb.toString();
    }
    
    public static ArrayList<Kmer> naiveExtendRight(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<String> terminators) {        
        ArrayList<Kmer> result = new ArrayList<>();
        
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
            
            if (best == null || terminators.contains(best.seq)) {
                break;
            }
            
            result.add(best);
            terminators.add(best.seq);
        }
        
        return result;
    }
    
    public static String naiveExtendLeft(String kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<String> terminators) {        
        StringBuilder sb = new StringBuilder(100);
        
        ArrayDeque<String> neighbors = graph.getPredecessors(kmer);
        String best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (String s : graph.getRightVariants(best)) {
                if (hasDepthRight(s, graph, maxTipLength)) {
                    sb.reverse();
                    return sb.toString();
                }
            }
            
            if (neighbors.size() == 1) {
                best = neighbors.peek();
            }
            else {
                best = null;
                for (String n : neighbors) {
                    if (hasDepthLeft(n, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            sb.reverse();
                            return sb.toString();
                        }
                    }
                }
            }
            
            if (best == null || terminators.contains(best)) {
                break;
            }
            
            sb.append(best.charAt(0));
            terminators.add(best);
        }
        
        sb.reverse();
        return sb.toString();
    }
    
    public static ArrayList<Kmer> naiveExtendLeft(Kmer kmer, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<String> terminators, boolean reverseResult) {        
        ArrayList<Kmer> result = new ArrayList<>();
        
        ArrayDeque<Kmer> neighbors = graph.getPredecessors(kmer);
        Kmer best = kmer;
        while (!neighbors.isEmpty()) {
            /** look for back branches*/
            for (Kmer s : graph.getRightVariants(best)) {
                if (hasDepthRight(s, graph, maxTipLength)) {
                    if (reverseResult) {
                        Collections.reverse(result);
                    }

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
                            if (reverseResult) {
                                Collections.reverse(result);
                            }

                            return result;
                        }
                    }
                }
            }
            
            if (best == null || terminators.contains(best.seq)) {
                break;
            }
            
            result.add(best);
            terminators.add(best.seq);
        }
        
        if (reverseResult) {
            Collections.reverse(result);
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
