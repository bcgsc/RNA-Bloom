/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
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
    
    public static Kmer greedyExtendRightOnce(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead) {
        LinkedList<Kmer> candidates = graph.getSuccessors(source);
        
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.peek();
            }
            else {
                Kmer cursor = candidates.pop();
                
                ArrayList<Kmer> path = new ArrayList<>(lookahead); 
                path.add(cursor);
                
                ArrayList<LinkedList<Kmer>> frontier = new ArrayList<>(lookahead);
                frontier.add(candidates);
                
                float bestCov = 0;
                ArrayList<Kmer> bestPath = path;
                
                while (!frontier.isEmpty()) {
                    if (path.size() < lookahead) {
                        candidates = graph.getSuccessors(cursor);
                        if (!candidates.isEmpty()) {
                            cursor = candidates.pop();
                            path.add(cursor);
                            frontier.add(candidates);
                            continue;
                        }
                    }
                    
                    float pathCov = getMedianKmerCoverage(path);
                    if (bestCov < pathCov) {
                        bestPath = new ArrayList<>(path);
                        bestCov = pathCov;
                    }

                    int i = path.size()-1;
                    while (i >= 0) {
                        candidates = frontier.get(i);
                        path.remove(i);
                        if (candidates.isEmpty()) {
                            frontier.remove(i);
                            --i;
                        }
                        else {
                            cursor = candidates.pop();
                            path.add(cursor);
                            break;
                        }
                    }
                }
                
                return bestPath.get(0);
            }
        }
    }
    
    public static Kmer greedyExtendLeftOnce(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead) {
        LinkedList<Kmer> candidates = graph.getPredecessors(source);
        
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.peek();
            }
            else {
                Kmer cursor = candidates.pop();
                
                ArrayList<Kmer> path = new ArrayList<>(lookahead); 
                path.add(cursor);
                
                ArrayList<LinkedList<Kmer>> frontier = new ArrayList<>(lookahead);
                frontier.add(candidates);
                
                float bestCov = 0;
                ArrayList<Kmer> bestPath = path;
                
                while (!frontier.isEmpty()) {
                    if (path.size() < lookahead) {
                        candidates = graph.getPredecessors(cursor);
                        if (!candidates.isEmpty()) {
                            cursor = candidates.pop();
                            path.add(cursor);
                            frontier.add(candidates);
                            continue;
                        }
                    }
                    
                    float pathCov = getMedianKmerCoverage(path);
                    if (bestCov < pathCov) {
                        bestPath = new ArrayList<>(path);
                        bestCov = pathCov;
                    }

                    int i = path.size()-1;
                    while (i >= 0) {
                        candidates = frontier.get(i);
                        path.remove(i);
                        if (candidates.isEmpty()) {
                            frontier.remove(i);
                            --i;
                        }
                        else {
                            cursor = candidates.pop();
                            path.add(cursor);
                            break;
                        }
                    }
                }
                
                return bestPath.get(0);
            }
        }
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
        LinkedList<Kmer> neighbors;
        for (int depth=0; depth < bound; ++depth) {
            neighbors = graph.getSuccessors(best);
            
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
                    /*right path intersects the left path */
                    String convergingKmer = best.seq;
                    ArrayList<Kmer> path = new ArrayList<>(bound);
                    for (Kmer kmer : leftPath) {
                        if (convergingKmer.equals(kmer.seq)) {
                            break;
                        }
                        else {
                            path.add(kmer);
                        }
                    }
                    Collections.reverse(rightPath);
                    path.addAll(rightPath);
                    return path;
                }
                else {
                    rightPath.add(best);
                }
            }
        }
        
        return null;
    }
    
    private static float getMedian(float[] arr) {
        float[] a = Arrays.copyOf(arr, arr.length);
        Arrays.sort(a);
        int halfLen = a.length/2;
        if (halfLen % 2 == 0) {
            return (a[halfLen-1] + a[halfLen])/2.0f;
        }
        
        return a[halfLen];
    }
    
    private static float getMedian(List<Float> arr) {
        ArrayList<Float> a = new ArrayList<>(arr);
        Collections.sort(a);
        int halfLen = a.size()/2;
        if (halfLen % 2 == 0) {
            return (a.get(halfLen-1) + a.get(halfLen))/2.0f;
        }
        
        return a.get(halfLen);
    }
    
    private static float rightGuidedMedianCoverageHelper(BloomFilterDeBruijnGraph graph, Kmer left, String guide) {
        int guideLen = guide.length();
        float[] covs = new float[guideLen+1];
        covs[guideLen] = left.count;
        
        String postfix = left.seq.substring(1);
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

    private static float leftGuidedMedianCoverageHelper(BloomFilterDeBruijnGraph graph, Kmer right, String guide) {
        int guideLen = guide.length();
        float[] covs = new float[guideLen+1];
        covs[0] = right.count;
        
        int kMinus1 = graph.getK()-1;
        String prefix = right.seq.substring(0,kMinus1);
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
    
    public static ArrayList<Kmer> correctMismatches(String seq, BloomFilterDeBruijnGraph graph, int lookahead, int mismatchesAllowed) {
        int seqLen = seq.length();
        int k = graph.getK();
        
        int numKmers = seqLen - k + 1;
        
        if (numKmers > 1) {
            char[] correctedSeq = seq.toCharArray();
            
            int mismatchesCorrected = 0;

            String currKmerSeq;
            Kmer prevKmer = graph.getKmer(seq.substring(0,k));
            String guide;
            Kmer bestKmer;
            float bestCov;
            float count;

            for (int i=k; i<seqLen-1; ++i) {
                currKmerSeq = graph.getSuffix(prevKmer.seq) + seq.charAt(i);
                guide = seq.substring(i+1, Math.min(i+1+lookahead, seqLen));

                bestKmer = null;
                bestCov = 0;

                for (Kmer s : graph.getSuccessors(prevKmer)) {
                    count = rightGuidedMedianCoverageHelper(graph, s, guide);
                    if (count > bestCov) {
                        bestKmer = s;
                        bestCov = count;
                    }
                }

                if (bestKmer == null) {
                    // undo all corrections
                    return graph.getKmers(seq);
                }
                else {
                    if (!currKmerSeq.equals(bestKmer.seq) && ++mismatchesCorrected > mismatchesAllowed) {
                        // too many mismatches, undo all corrections
                        return graph.getKmers(seq);
                    }

                    correctedSeq[i] = graph.getLastBase(bestKmer.seq);
                    prevKmer = bestKmer;
                }
            }
            
            bestKmer = greedyExtendRightOnce(graph, prevKmer, lookahead);
            if (bestKmer != null) {
                correctedSeq[seqLen-1] = graph.getLastBase(bestKmer.seq);
            }
            
            /** correct mismatches in first kmer of the sequence*/
            
            /** Get the k-th kmer or the last kmer, whichever is more to the left*/
            int i = Math.min(2*k, seqLen);
            prevKmer = graph.getKmer(seq.substring(i-k,i));
            
            for (i=i-k-1; i>0; --i) {
                currKmerSeq = seq.charAt(i) + graph.getPrefix(prevKmer.seq);
                guide = seq.substring(Math.max(0, i-lookahead), i);
                
                bestKmer = null;
                bestCov = 0;

                for (Kmer s : graph.getPredecessors(prevKmer)) {
                    count = leftGuidedMedianCoverageHelper(graph, s, guide);
                    if (count > bestCov) {
                        bestKmer = s;
                        bestCov = count;
                    }
                }

                if (bestKmer == null) {
                    // undo all corrections
                    return graph.getKmers(seq);
                }
                else {
                    if (!currKmerSeq.equals(bestKmer.seq) && ++mismatchesCorrected > mismatchesAllowed) {
                        // too many mismatches, undo all corrections
                        return graph.getKmers(seq);
                    }

                    correctedSeq[i] = graph.getFirstBase(bestKmer.seq);
                    prevKmer = bestKmer;
                }
            }
            
            bestKmer = greedyExtendLeftOnce(graph, prevKmer, lookahead);
            if (bestKmer != null) {
                correctedSeq[0] = graph.getFirstBase(bestKmer.seq);
            }
            
            return graph.getKmers(new String(correctedSeq));
        }
        else {
            return graph.getKmers(seq);
        }
    }

    public static String assemble(ArrayList<Kmer> kmers) {
        
        String first = kmers.get(0).seq;
        int k = first.length();
        int lastIndex = k - 1;
        
        StringBuilder sb = new StringBuilder(k + kmers.size() - 1);
        sb.append(first.substring(0, lastIndex));
        
        for (Kmer kmer : kmers) {
            sb.append(kmer.seq.charAt(lastIndex));
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
                
        LinkedList<Float> window = new LinkedList<>();
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
        ArrayList<Kmer> path = null;
        
        Kmer best = seed;
        path = greedyExtend(best, graph, lookahead);
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
    
    public static String assembleFragment(String leftRead, String rightRead, BloomFilterDeBruijnGraph graph, int mismatchesAllowed, int bound, int lookahead, int minOverlap) {
        
        // overlap before finding path
        String overlapped = overlapMaximally(leftRead, rightRead, minOverlap);
        if (overlapped != null) {
            /**@TODO Check whether overlap is a valid path in DBG*/
            return assemble(correctMismatches(overlapped, graph, lookahead, mismatchesAllowed));
        }
        
        ArrayList<Kmer> leftKmers = correctMismatches(leftRead, graph, lookahead, mismatchesAllowed);
        ArrayList<Kmer> rightKmers = correctMismatches(rightRead, graph, lookahead, mismatchesAllowed);
        
        ArrayList<Kmer> pathKmers = getMaxCoveragePath(graph, leftKmers.get(leftKmers.size()-1), rightKmers.get(0), bound, lookahead);
        
        if (pathKmers == null || pathKmers.isEmpty()) {
            return "";
        }
        
        return assembleFirstBase(leftKmers) + assemble(pathKmers) + assembleLastBase(rightKmers);
    }
    
    public static boolean hasDepthRight(Kmer source, BloomFilterDeBruijnGraph graph, int depth) {
        LinkedList<LinkedList> frontier = new LinkedList<>();
        LinkedList<Kmer> alts = graph.getSuccessors(source);
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
        LinkedList<LinkedList> frontier = new LinkedList<>();
        LinkedList<Kmer> alts = graph.getPredecessors(source);
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
    
    public static String naiveExtendRight(Kmer source, BloomFilterDeBruijnGraph graph, int maxTipLength) {        
        StringBuilder sb = new StringBuilder(50);
        int i = graph.getK() - 1;
        
        LinkedList<Kmer> neighbors = graph.getSuccessors(source);
        while (!neighbors.isEmpty()) {
            if (neighbors.size() == 1) {
                sb.append(neighbors.peek().seq.charAt(i));
            }
            else {
                Kmer best = null;
                for (Kmer n : neighbors) {
                    if (hasDepthRight(n, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            best = null;
                            break;
                        }
                    }
                }
                
                if (best == null) {
                    break;
                }
                
                sb.append(best.seq.charAt(i));
            }
        }
        
        return sb.toString();
    }
    
    public static String naiveExtendLeft(Kmer source, BloomFilterDeBruijnGraph graph, int maxTipLength) {        
        StringBuilder sb = new StringBuilder(50);
        int i = graph.getK() - 1;
        
        LinkedList<Kmer> neighbors = graph.getPredecessors(source);
        while (!neighbors.isEmpty()) {
            if (neighbors.size() == 1) {
                sb.append(neighbors.peek().seq.charAt(i));
            }
            else {
                Kmer best = null;
                for (Kmer n : neighbors) {
                    if (hasDepthLeft(n, graph, maxTipLength)) {
                        if (best == null) {
                            best = n;
                        }
                        else {
                            // too many good branches
                            best = null;
                            break;
                        }
                    }
                }
                
                if (best == null) {
                    break;
                }
                
                sb.append(best.seq.charAt(i));
            }
        }
        
        return sb.toString();
    }
}
