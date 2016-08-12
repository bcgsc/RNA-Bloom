/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;

/**
 *
 * @author gengar
 */
public final class GraphUtils {
    
    public static Kmer greedyExtendRightOnce(BloomFilterDeBruijnGraph graph, Kmer source, int lookahead) {
        ArrayList<Kmer> candidates = graph.getSuccessors(source);
        
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.get(0);
            }
            else {
                ArrayList<Kmer> alts = graph.getSuccessors(source);
                Kmer cursor = alts.remove(alts.size()-1);
                
                ArrayList<Kmer> path = new ArrayList<>(lookahead); 
                path.add(cursor);
                
                ArrayList<ArrayList<Kmer>> frontier = new ArrayList<>(lookahead);
                frontier.add(alts);
                
                float bestCov = 0;
                int bestLen = 1;
                ArrayList<Kmer> bestPath = path;
                
                while (!frontier.isEmpty()) {
                    if (path.size() < lookahead) {
                        alts = graph.getSuccessors(cursor);
                        if (!alts.isEmpty()) {
                            cursor = alts.remove(alts.size()-1);
                            path.add(cursor);
                            frontier.add(alts);
                            continue;
                        }
                    }
                    
                    float pathCov = graph.getMedianKmerCoverage(path);
                    int pathLen = path.size();
                    if (bestLen < pathLen || bestCov < pathCov) {
                        bestPath = new ArrayList<>(path);
                        bestCov = pathCov;
                        bestLen = pathLen;
                    }

                    int i = path.size()-1;
                    while (i >= 0) {
                        alts = frontier.get(i);
                        path.remove(i);
                        if (alts.isEmpty()) {
                            frontier.remove(i);
                            --i;
                        }
                        else {
                            cursor = alts.remove(alts.size()-1);
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
        ArrayList<Kmer> candidates = graph.getPredecessors(source);
        
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.get(0);
            }
            else {
                ArrayList<Kmer> alts = graph.getPredecessors(source);
                Kmer cursor = alts.remove(alts.size()-1);
                
                ArrayList<Kmer> path = new ArrayList<>(lookahead); 
                path.add(cursor);
                
                ArrayList<ArrayList<Kmer>> frontier = new ArrayList<>(lookahead);
                frontier.add(alts);
                
                float bestCov = 0;
                int bestLen = 1;
                ArrayList<Kmer> bestPath = path;
                
                while (!frontier.isEmpty()) {
                    if (path.size() < lookahead) {
                        alts = graph.getPredecessors(cursor);
                        if (!alts.isEmpty()) {
                            cursor = alts.remove(alts.size()-1);
                            path.add(cursor);
                            frontier.add(alts);
                            continue;
                        }
                    }
                    
                    float pathCov = graph.getMedianKmerCoverage(path);
                    int pathLen = path.size();
                    if (bestLen < pathLen || bestCov < pathCov) {
                        bestPath = new ArrayList<>(path);
                        bestCov = pathCov;
                        bestLen = pathLen;
                    }

                    int i = path.size()-1;
                    while (i >= 0) {
                        alts = frontier.get(i);
                        path.remove(i);
                        if (alts.isEmpty()) {
                            frontier.remove(i);
                            --i;
                        }
                        else {
                            cursor = alts.remove(alts.size()-1);
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
        ArrayList<Kmer> neighbors;
        for (int depth=0; depth < bound; ++depth) {
            neighbors = graph.getSuccessors(best);
            
            if (neighbors.isEmpty()) {
                break;
            }
            else {
                if (neighbors.size() == 1) {
                    best = neighbors.get(0);
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
                    best = neighbors.get(0);
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
    
    private static float getMedianCoverageHelper(BloomFilterDeBruijnGraph graph, Kmer left, String guide) {
        ArrayList<Kmer> kmers = new ArrayList<>(guide.length()+1);
        kmers.add(left);
        
        String prevPostfix = left.seq.substring(1);
        int guideLen = guide.length();
        String kmer;
        float count;
        for (int i=0; i<guideLen; ++i) {
            kmer = prevPostfix + guide.charAt(i);
            count = graph.getCount(kmer);
            if (count > 0) {
                kmers.add(new Kmer(kmer, count));
                prevPostfix = kmer.substring(1);
            }
            else {
                // not a valid sequence
                return 0;
            }
        }
        
        return graph.getMedianKmerCoverage(kmers);
    }
    
    public static ArrayList<Kmer> correctMismatches(String seq, BloomFilterDeBruijnGraph graph, int lookahead, int mismatchesAllowed) {
        int seqLen = seq.length();
        int k = graph.getK();
        
        int numKmers = seqLen - k + 1;        
        ArrayList<Kmer> correctedKmers = new ArrayList<>(numKmers);
        
        int mismatchesCorrected = 0;
        
        String currKmerSeq;
        String prevKmerSeq = seq.substring(0,k);
        Kmer prevKmer = new Kmer(prevKmerSeq, graph.getCount(prevKmerSeq));
        String guide;
        Kmer bestKmer;
        float bestCov;
        float count;
                
        for (int i=0; i<numKmers-1; ++i) {
            currKmerSeq = prevKmer.seq.substring(1) + seq.charAt(i+k);
            guide = seq.substring(i+k+1, Math.min(i+k+1+lookahead, seqLen));
            
            bestKmer = null;
            bestCov = 0;
            
            for (Kmer s : graph.getSuccessors(prevKmer)) {
                count = getMedianCoverageHelper(graph, s, guide);
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
                if (! currKmerSeq.equals(bestKmer.seq)) {
                    ++mismatchesCorrected;
                    if (mismatchesCorrected > mismatchesAllowed) {
                        return graph.getKmers(seq);
                    }
                }
                
                correctedKmers.add(bestKmer);
                prevKmer = bestKmer;
            }
        }
        
        /**@TODO correct mismatches in first kmer*/
        
        return correctedKmers;
    }
}
