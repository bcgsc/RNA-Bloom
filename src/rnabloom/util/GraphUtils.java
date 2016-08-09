/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.ArrayList;
import java.util.HashMap;
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
                    
                    float pathCov = graph.getMeanKmerCoverage(path);
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
                    
                    float pathCov = graph.getMeanKmerCoverage(path);
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
        ArrayList<Kmer> path = new ArrayList<>(bound);
        
        /* extend right */
        Kmer best = left;
        ArrayList<Kmer> successors;
        for (int depth=0; depth < bound; ++depth) {
            successors = graph.getSuccessors(best);
            
            if (successors.isEmpty()) {
                break;
            }
            else {
                if (successors.size() == 1) {
                    best = successors.get(0);
                }
                else {
                    /**@TODO */
                }
            }
        }
        
        return path;
    }
    
}
