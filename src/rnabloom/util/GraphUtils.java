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
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;
import static rnabloom.util.SeqUtils.kmerize;
import static rnabloom.util.SeqUtils.kmerizeArrayList;
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
    
    private static float rightGuidedMedianCoverage(BloomFilterDeBruijnGraph graph, String source, String guide) {
        int guideLen = guide.length();
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
    
    private static float leftGuidedMedianCoverage(BloomFilterDeBruijnGraph graph, String source, String guide) {
        int guideLen = guide.length();
        if (guideLen == 0) {
            return graph.getCount(source);
        }
        
        float[] covs = new float[guideLen+1];
        covs[0] = graph.getCount(source);
        
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
    
    public static String correctMismatches(String seq, BloomFilterDeBruijnGraph graph, int lookahead, int mismatchesAllowed) {
        
        StringBuilder sb = new StringBuilder(seq);
        int seqLen = seq.length();
        int k = graph.getK();
        int numKmers = seqLen-k+1;
        
        // correct from start
        for (int i=0; i<numKmers; ++i) {
            int end = i+k;
            String kmer = sb.substring(i, end);
            LinkedList<String> variants = graph.getRightVariants(kmer);
            if (!variants.isEmpty()) {
                String guide = sb.substring(end, Math.min(end+lookahead, seqLen));
                float bestCov = rightGuidedMedianCoverage(graph, kmer, guide);
                
                boolean corrected = false;
                for (String v : variants) {
                    float cov = rightGuidedMedianCoverage(graph, v, guide);
                    if (cov > bestCov) {
                        bestCov = cov;
                        sb.setCharAt(end-1, v.charAt(k-1));
                        corrected = true;
                    }
                }
                
                if (corrected) {
                    --mismatchesAllowed;
                    if (mismatchesAllowed < 0) {
                        // too many mismatches
                        return seq;
                    }
                }
            }
        }
        
        // correct from end
        for (int i=seqLen-k; i>=0; --i) {
            String kmer = sb.substring(i, i+k);
            LinkedList<String> variants = graph.getLeftVariants(kmer);
            if (!variants.isEmpty()) {
                String guide = sb.substring(Math.max(0, i-lookahead), i);
                float bestCov = leftGuidedMedianCoverage(graph, kmer, guide);
                
                boolean corrected = false;
                for (String v : variants) {
                    float cov = leftGuidedMedianCoverage(graph, v, guide);
                    if (cov > bestCov) {
                        bestCov = cov;
                        sb.setCharAt(i, v.charAt(0));
                        corrected = true;
                    }
                }
                
                if (corrected) {
                    --mismatchesAllowed;
                    if (mismatchesAllowed < 0) {
                        // too many mismatches
                        return seq;
                    }
                }
            }
        }
        
        return sb.toString();
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
    
    public static String assembleFragment(String leftRead, String rightRead, BloomFilterDeBruijnGraph graph, int mismatchesAllowed, int bound, int lookahead, int minOverlap) {
        
        // overlap before finding path
        String overlapped = overlapMaximally(leftRead, rightRead, minOverlap);
        if (overlapped != null) {
            /**@TODO Check whether overlap is a valid path in DBG*/
            return correctMismatches(overlapped, graph, lookahead, mismatchesAllowed);
        }
        
        leftRead = correctMismatches(leftRead, graph, lookahead, mismatchesAllowed);
        rightRead = correctMismatches(rightRead, graph, lookahead, mismatchesAllowed);
        
        int k = graph.getK();
        int leftReadLength = leftRead.length();
        Kmer leftReadRightmostKmer = graph.getKmer(leftRead.substring(leftReadLength-k, leftReadLength));
        Kmer rightReadLeftmostKmer = graph.getKmer(rightRead.substring(0, k));
        
        ArrayList<Kmer> pathKmers = getMaxCoveragePath(graph, leftReadRightmostKmer, rightReadLeftmostKmer, bound, lookahead);
        
        if (pathKmers == null || pathKmers.isEmpty()) {
            return "";
        }
        
        String leftWing, rightWing;
        if (leftReadLength == k) {
            leftWing = leftRead.substring(0, 1);
        }
        else {
            leftWing = leftRead.substring(0, leftReadLength-k+1);
        }
        
        int rightReadLength = leftRead.length();
        if (rightReadLength == k) {
            rightWing = rightRead.substring(k-1, k);
        }
        else {
            rightWing = rightRead.substring(k-1);
        }
        
        return leftWing + assemble(pathKmers) + rightWing;
    }
    
    public static String extendWithPairedKmers(String fragment, BloomFilterDeBruijnGraph graph, int lookAhead) {
        int distance = graph.getPairedKmerDistance();
        ArrayList<String> kmers = kmerizeArrayList(fragment, graph.getK());
        HashSet<String> fragmentKmers = new HashSet<>(kmers);
        
        /** extend right*/
        String best = kmers.get(kmers.size()-1);
        LinkedList<String> neighbors = graph.getSuccessors(best);
        while (!neighbors.isEmpty()) {
            best = null;
            if (neighbors.size() == 1) {
                best = neighbors.peek();
            }
            else {
                // >1 neighbors
                
                LinkedList<String> fragmentNeighbors = new LinkedList<>();
                for (String n : neighbors) {
                    if (graph.lookupFragmentKmer(n)) {
                        fragmentNeighbors.add(n);
                    }
                }
                
                if (fragmentNeighbors.isEmpty()) {
                    break;
                }
                else if (fragmentNeighbors.size() == 1) {
                    best = fragmentNeighbors.peek();
                }
                else {
                    // >1 fragment neighbors
                    
                    LinkedList<String> pairedNeighbors = new LinkedList<>();
                    
                    if (kmers.size() >= distance) {
                        String partner = kmers.get(kmers.size()-distance);
                        for (String n : fragmentNeighbors) {
                            if (graph.lookupPairedKmers(partner, n)) {
                                pairedNeighbors.add(n);
                            }
                        }
                    }
                    
                    if (pairedNeighbors.isEmpty()) {
                        // check depth
                        LinkedList<String> longFragmentNeighbors = new LinkedList<>();
                        for (String n : fragmentNeighbors) {
                            if (hasFragmentDepthRight(n, graph, lookAhead)) {
                                longFragmentNeighbors.add(n);
                            }
                        }
                        
                        if (longFragmentNeighbors.isEmpty()) {
                            break;
                        }
                        else if (longFragmentNeighbors.size() == 1) {
                            best = longFragmentNeighbors.peek();
                        }
                        else {
                            // >1 long fragment neighbors
                            break;
                        }
                    }
                    else if (pairedNeighbors.size() == 1) {
                        best = pairedNeighbors.peek();
                    }
                    else {
                        // >1 paired neighbors
                        break;
                    }
                }
            }
            
            if (best == null || fragmentKmers.contains(best)) {
                break;
            }
            
            kmers.add(best);
            fragmentKmers.add(best);
            neighbors = graph.getSuccessors(best);
        }
        
        Collections.reverse(kmers);
        
        /** extend left*/
        best = kmers.get(kmers.size()-1);
        neighbors = graph.getPredecessors(best);
        while (!neighbors.isEmpty()) {
            best = null;
            if (neighbors.size() == 1) {
                best = neighbors.peek();
            }
            else {
                // >1 neighbors
                
                LinkedList<String> fragmentNeighbors = new LinkedList<>();
                for (String n : neighbors) {
                    if (graph.lookupFragmentKmer(n)) {
                        fragmentNeighbors.add(n);
                    }
                }
                
                if (fragmentNeighbors.isEmpty()) {
                    break;
                }
                else if (fragmentNeighbors.size() == 1) {
                    best = fragmentNeighbors.peek();
                }
                else {
                    // >1 fragment neighbors
                    
                    LinkedList<String> pairedNeighbors = new LinkedList<>();
                    if (kmers.size() >= distance) {
                        String partner = kmers.get(kmers.size()-distance);
                        for (String n : fragmentNeighbors) {
                            if (graph.lookupPairedKmers(n, partner)) {
                                pairedNeighbors.add(n);
                            }
                        }
                    }
                    
                    if (pairedNeighbors.isEmpty()) {
                        // check depth
                        LinkedList<String> longFragmentNeighbors = new LinkedList<>();
                        for (String n : fragmentNeighbors) {
                            if (hasFragmentDepthLeft(n, graph, lookAhead)) {
                                longFragmentNeighbors.add(n);
                            }
                        }
                        
                        if (longFragmentNeighbors.isEmpty()) {
                            break;
                        }
                        else if (longFragmentNeighbors.size() == 1) {
                            best = longFragmentNeighbors.peek();
                        }
                        else {
                            // >1 long fragment neighbors
                            break;
                        }
                    }
                    else if (pairedNeighbors.size() == 1) {
                        best = pairedNeighbors.peek();
                    }
                    else {
                        // >1 paired neighbors
                        break;
                    }
                }
            }
            
            if (best == null || fragmentKmers.contains(best)) {
                break;
            }
            
            kmers.add(best);
            fragmentKmers.add(best);
            neighbors = graph.getPredecessors(best);
        }
        
        Collections.reverse(kmers);
        
        return assembleString(kmers);
    }
    
    private static boolean hasFragmentDepthRight(String source, BloomFilterDeBruijnGraph graph, int depth) {
        LinkedList<LinkedList> frontier = new LinkedList<>();
        LinkedList<String> alts = new LinkedList<>();
        for (String s : graph.getSuccessors(source)) {
            if (graph.lookupFragmentKmer(s)) {
                alts.add(s);
            }
        }
        frontier.add(alts);
        
        while (!frontier.isEmpty()) {
            alts = frontier.peekLast();
            if (alts.isEmpty()) {
                frontier.removeLast();
            }
            else {
                String a = alts.pop();
                alts = new LinkedList<>();
                for (String s : graph.getSuccessors(a)) {
                    if (graph.lookupFragmentKmer(s)) {
                        alts.add(s);
                    }
                }
                frontier.add(alts);
            }

            if (frontier.size() >= depth) {
                return true;
            }
        }
        
        return false;
    }
    
    private static boolean hasFragmentDepthLeft(String source, BloomFilterDeBruijnGraph graph, int depth) {
        LinkedList<LinkedList> frontier = new LinkedList<>();
        LinkedList<String> alts = new LinkedList<>();
        for (String s : graph.getPredecessors(source)) {
            if (graph.lookupFragmentKmer(s)) {
                alts.add(s);
            }
        }
        frontier.add(alts);
        
        while (!frontier.isEmpty()) {
            alts = frontier.peekLast();
            if (alts.isEmpty()) {
                frontier.removeLast();
            }
            else {
                String a = alts.pop();
                alts = new LinkedList<>();
                for (String s : graph.getPredecessors(a)) {
                    if (graph.lookupFragmentKmer(s)) {
                        alts.add(s);
                    }
                }
                frontier.add(alts);
            }

            if (frontier.size() >= depth) {
                return true;
            }
        }
        
        return false;
    }
    
    private static boolean hasDepthRight(Kmer source, BloomFilterDeBruijnGraph graph, int depth) {
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
    
    private static boolean hasDepthLeft(Kmer source, BloomFilterDeBruijnGraph graph, int depth) {
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
    
    public static String naiveExtend(String fragment, BloomFilterDeBruijnGraph graph, int maxTipLength) {
        String[] kmers = kmerize(fragment, graph.getK());
        int numKmers = kmers.length;
        HashSet<String> fragmentKmers = new HashSet<>(2*numKmers);
        String extended = "";
        
        ArrayList<Kmer> extension = naiveExtendLeft(graph.getKmer(kmers[0]), graph, maxTipLength, fragmentKmers);
        if (!extension.isEmpty()) {
            extended = assembleFirstBase(extension);
        }
        
        extended += fragment;
        
        extension = naiveExtendRight(graph.getKmer(kmers[numKmers-1]), graph, maxTipLength, fragmentKmers);
        if (!extension.isEmpty()) {
            extended += assembleLastBase(extension);
        }
        
        return extended;
    }
    
    
    private static ArrayList<Kmer> naiveExtendRight(Kmer source, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<String> terminators) {        
        ArrayList<Kmer> extension = new ArrayList<>(100);
        
        LinkedList<Kmer> neighbors = graph.getSuccessors(source);
        Kmer best;
        while (!neighbors.isEmpty()) {
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
                            best = null;
                            break;
                        }
                    }
                }
                
                /**@TODO look for back branches*/
                
            }
            
            if (best == null || terminators.contains(best.seq)) {
                break;
            }

            extension.add(best);
            terminators.add(best.seq);
        }
        
        return extension;
    }
    
    private static ArrayList<Kmer> naiveExtendLeft(Kmer source, BloomFilterDeBruijnGraph graph, int maxTipLength, HashSet<String> terminators) {        
        ArrayList<Kmer> extension = new ArrayList<>(100);
        
        LinkedList<Kmer> neighbors = graph.getPredecessors(source);
        Kmer best;
        while (!neighbors.isEmpty()) {
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
                            best = null;
                            break;
                        }
                    }
                }
                
                /**@TODO look for back branches*/
                
            }
            
            if (best == null || terminators.contains(best.seq)) {
                break;
            }

            extension.add(best);
            terminators.add(best.seq);
        }
        
        Collections.reverse(extension);
        
        return extension;
    }
}
