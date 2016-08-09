/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import rnabloom.graph.BloomFilterDeBruijnGraph;
import rnabloom.graph.BloomFilterDeBruijnGraph.Kmer;

/**
 *
 * @author gengar
 */
public final class GraphUtils {
    
    public static Kmer greedyExtendRightOnce(BloomFilterDeBruijnGraph graph, String source, int lookahead) {
        ArrayList<Kmer> candidates = graph.getSuccessors(source);
        
        if (candidates.isEmpty()) {
            return null;
        }
        else {
            if (candidates.size() == 1) {
                return candidates.get(0);
            }
            else {
                Kmer best = null;
                
                /**@TODO */
                
                return best;
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
