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
package rnabloom.olc;

import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.TransitiveReduction;
import org.jgrapht.alg.connectivity.KosarajuStrongConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import rnabloom.io.ExtendedPafRecord;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.PafReader;
import rnabloom.io.PafRecord;
import static rnabloom.util.SeqUtils.stringToBytes;
import static rnabloom.util.SeqUtils.bytesToString;
import static rnabloom.util.SeqUtils.reverseComplement;


/**
 *
 * @author kmnip
 */
public class Layout {
    
    private final static Pattern CIGAR_OP_PATTERN = Pattern.compile("(\\d+)([MIDNSHPX=])");
    private DefaultDirectedGraph<String, OverlapEdge> graph;
    private String overlapPafPath;
    private String seqFastaPath;
    private boolean stranded;
    private int maxEdgeClip = 100;
    private float minAlnId = 0.50f;
    private int maxIndelSize = 20;
    private int minOverlapMatches = 200;
    private boolean cutRevCompArtifact = false;
    private int minNumAltReads = 0;
    
    public Layout(String seqFile, String pafFile, boolean stranded, int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact, int minSeqDepth) {
        this.graph = new DefaultDirectedGraph<>(OverlapEdge.class);
        this.overlapPafPath = pafFile;
        this.seqFastaPath = seqFile;
        this.stranded = stranded;
        this.maxEdgeClip = maxEdgeClip;
        this.minAlnId = minAlnId;
        this.minOverlapMatches = minOverlapMatches;
        this.maxIndelSize = maxIndelSize;
        this.cutRevCompArtifact = cutRevCompArtifact;
        this.minNumAltReads = minSeqDepth - 1;
    }
    
    private class OverlapEdge extends DefaultEdge implements Comparable<OverlapEdge> {
        public int sourceStart, sourceEnd, sinkStart, sinkEnd;
        
        public OverlapEdge(int sourceStart, int sourceEnd, int sinkStart, int sinkEnd) {
            this.sourceStart = sourceStart;
            this.sourceEnd = sourceEnd;
            this.sinkStart = sinkStart;
            this.sinkEnd = sinkEnd;
        }

        @Override
        public int compareTo(OverlapEdge o) {
            // an edge with a larger overlap is "less"
            return Math.max(o.sinkEnd-o.sinkStart, o.sourceEnd-o.sourceStart) - Math.max(sinkEnd-sinkStart, sourceEnd-sourceStart);
        }
    }

    private boolean hasGoodOverlap(PafRecord r) {
        return r.numMatch / (float)(r.qEnd - r.qStart) >= minAlnId &&
               r.numMatch / (float)(r.tEnd - r.tStart) >= minAlnId;
    }
    
    private boolean hasAlignment(ExtendedPafRecord record) {
        return record.cigar != null && record.nm >= 0;
    }
    
    private boolean hasGoodAlignment(ExtendedPafRecord record) {
        int numMatch = 0;
        int numDel = 0;
        int numIns = 0;
                
        Matcher m = CIGAR_OP_PATTERN.matcher(record.cigar);
        while (m.find()) {
            int opSize = Integer.parseInt(m.group(1));
            char op = m.group(2).charAt(0);
            switch (op) {
                case 'M':
                    numMatch += opSize;
                    break;
                case 'I':
                    if (opSize > maxIndelSize) {
                        return false;
                    }
                    numIns += opSize;
                    break;
                case 'D':
                    if (opSize > maxIndelSize) {
                        return false;
                    }
                    numDel += opSize;
                    break;
            }
        }
        
        float alnId = (numMatch - record.nm)/(float)(numMatch + numDel + numIns);
        
        return alnId >= minAlnId;
    }
    
    private boolean hasReverseComplementArtifact(ExtendedPafRecord r) {
        if (r.qName.equals(r.tName) && r.reverseComplemented) {
            if (r.qStart <= maxEdgeClip || r.qLen - r.qEnd <= maxEdgeClip) {
                return true;
            }
        }
        
        return false;
    }
    
    private int getReverseComplementArtifactCutIndex(ExtendedPafRecord r) {
        if (r.qName.equals(r.tName) && r.reverseComplemented) {
            if (r.qStart <= maxEdgeClip) {
                // left edge has reverse complement
                if (r.qStart == r.tStart && r.qEnd == r.tEnd) {
                    // palindrome
                    return r.qStart + (r.qEnd - r.qStart)/2;
                }
                else {
                    return r.qEnd -1;
                }
            }
            else if (r.qLen - r.qEnd <= maxEdgeClip) {
                // right edge has reverse complement
                if (r.qStart == r.tStart && r.qEnd == r.tEnd) {
                    // palindrome
                    return r.qStart + (r.qEnd - r.qStart)/2;
                }
                else {
                    return r.tStart -1;
                }
            }
        }
        
        return -1;
    }
    
    private boolean isContainmentPafRecord(ExtendedPafRecord r) {
        return ((r.qStart <= maxEdgeClip && r.qLen - r.qEnd <= maxEdgeClip) ||
            (r.tStart <= maxEdgeClip && r.tLen - r.tEnd <= maxEdgeClip));
    }
    
    private boolean isStrandedContainmentPafRecord(ExtendedPafRecord r) {
        return !r.reverseComplemented && isContainmentPafRecord(r);
    }
        
    private boolean isDovetailPafRecord(ExtendedPafRecord r) {
        if (r.numMatch >= minOverlapMatches) {

            if (r.reverseComplemented) {
                return (r.qEnd >= r.qLen - maxEdgeClip && r.tEnd >= r.tLen - maxEdgeClip && r.qStart > r.tLen - r.tEnd) ||
                        (r.tStart <= maxEdgeClip && r.qStart <= maxEdgeClip && r.qLen - r.qStart > r.tStart);
            }
            else {
                return (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip && r.qStart > r.tStart) ||
                        (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip && r.tStart > r.qStart);
            }

        }
        
        return false;
    }
    
    private boolean isStrandedDovetailPafRecord(ExtendedPafRecord r) {
        if (!r.reverseComplemented) {
            if (r.numMatch >= minOverlapMatches &&
                ((r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) ||
                    (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip))) {
                return true;
            }
        }
        
        return false;
    }
    
    private static String getVertexName(String vid) {
        return vid.substring(0, vid.length()-1);
    }
    
    private static char geVertexSign(String vid) {
        return vid.charAt(vid.length()-1);
    }
    
    private static boolean isVertexSignReverseComplement(String vid) {
        return geVertexSign(vid) == '-';
    }
    
    private static char getReverseComplementSign(char sign) {
        switch (sign) {
            case '+':
                return '-';
            case '-':
                return '+';
        }
        
        return '0';
    }
    
    private static String getReverseComplementID(String vid) {
        int lastIndex = vid.length()-1;
        return vid.substring(0, lastIndex) + getReverseComplementSign(vid.charAt(lastIndex));
    }
    
    private OverlapEdge getReverseComplementEdge(OverlapEdge e) {
        String source = graph.getEdgeSource(e);
        String sink = graph.getEdgeTarget(e);
        return graph.getEdge(getReverseComplementID(sink), getReverseComplementID(source));
    }
    
    private int reduceTransitively() {
        int numEdgesRemoved = 0;
        
        KosarajuStrongConnectivityInspector<String, OverlapEdge> ci = new KosarajuStrongConnectivityInspector<>(graph);
        
        // Perform transitive reduction on each biconnected component; 
        // this routine should use less memory than reducing the entire graph.
        for (Graph<String, OverlapEdge> cc : ci.getStronglyConnectedComponents()) {
            Set<OverlapEdge> edges = cc.edgeSet();
            if (!edges.isEmpty()) {
                edges = new HashSet<>(edges);

                // create directed graph for the connected component
                DefaultDirectedGraph<String, OverlapEdge> g = new DefaultDirectedGraph<>(OverlapEdge.class);
                for (OverlapEdge e : edges) {
                    String source = graph.getEdgeSource(e);
                    String target = graph.getEdgeTarget(e);
                    g.addVertex(source);
                    g.addVertex(target);
                    g.addEdge(source, target, e);
                }

                TransitiveReduction.INSTANCE.reduce(g);

                Set<OverlapEdge> reducedGraphEdges = g.edgeSet();

                if (edges.size() > reducedGraphEdges.size()) {
                    edges.removeAll(reducedGraphEdges);
                    numEdgesRemoved += edges.size();
                    graph.removeAllEdges(edges);
                }
            }
        }
        
        return numEdgesRemoved;
    }
        
    private void resolveJunctions() {
        ArrayList<OverlapEdge> edges = new ArrayList<>(graph.edgeSet());
        Collections.sort(edges);
        for (OverlapEdge bestEdge : edges) {
            if (graph.containsEdge(bestEdge)) {
                String source = graph.getEdgeSource(bestEdge);
                ArrayDeque<OverlapEdge> edgesToRemove = new ArrayDeque<>();
                for (OverlapEdge e : graph.outgoingEdgesOf(source)) {
                    if (!e.equals(bestEdge)) {
                        edgesToRemove.add(e);
                    }
                }
                graph.removeAllEdges(edgesToRemove);
                
                String target = graph.getEdgeTarget(bestEdge);
                edgesToRemove = new ArrayDeque<>();
                for (OverlapEdge e : graph.incomingEdgesOf(target)) {
                    if (!e.equals(bestEdge)) {
                        edgesToRemove.add(e);
                    }
                }
                graph.removeAllEdges(edgesToRemove);
            }
        }
    }
    
    private void resolveJunctions(HashSet<String> names, boolean strandSpecific) {
        for (String name : names) {
            name += '+';
            
            if (graph.outDegreeOf(name) > 1) {
                Set<OverlapEdge> edges = graph.outgoingEdgesOf(name);
                ArrayDeque<OverlapEdge> edgesToRemove = new ArrayDeque<>(edges.size()-1);
                
                // find the best edge
                Iterator<OverlapEdge> itr = edges.iterator();
                OverlapEdge bestEdge = itr.next();
                int bestOverlap = bestEdge.sourceEnd - bestEdge.sourceStart;
                while (itr.hasNext()) {
                    OverlapEdge e = itr.next();
                    int overlap = e.sourceEnd - e.sourceStart;
                    if (overlap > bestOverlap) {
                        edgesToRemove.add(bestEdge);
                        
                        bestOverlap = overlap;
                        bestEdge = e;
                    }
                    else {
                        edgesToRemove.add(e);
                    }
                }
                
                if (!strandSpecific) {
                    ArrayDeque<OverlapEdge> reEdgesToRemove = new ArrayDeque<>(edgesToRemove.size());
                    for (OverlapEdge e : edgesToRemove) {
                        OverlapEdge re = getReverseComplementEdge(e);
                        if (re != null) {
                            reEdgesToRemove.add(re);
                        }                        
                    }
                    graph.removeAllEdges(reEdgesToRemove);
                }
                                
                // remove other edges from graph
                graph.removeAllEdges(edgesToRemove);
            }
            
            if (graph.inDegreeOf(name) > 1) {
                Set<OverlapEdge> edges = graph.incomingEdgesOf(name);
                ArrayDeque<OverlapEdge> edgesToRemove = new ArrayDeque<>(edges.size()-1);
                
                // find the best edge
                Iterator<OverlapEdge> itr = edges.iterator();
                OverlapEdge bestEdge = itr.next();
                int bestOverlap = bestEdge.sinkEnd - bestEdge.sinkStart;
                while (itr.hasNext()) {
                    OverlapEdge e = itr.next();
                    int overlap = e.sinkEnd - e.sinkStart;
                    if (overlap > bestOverlap) {
                        edgesToRemove.add(bestEdge);
                        
                        bestOverlap = overlap;
                        bestEdge = e;
                    }
                    else {
                        edgesToRemove.add(e);
                    }
                }
                
                if (!strandSpecific) {
                    ArrayDeque<OverlapEdge> reEdgesToRemove = new ArrayDeque<>(edgesToRemove.size());
                    for (OverlapEdge e : edgesToRemove) {
                        OverlapEdge re = getReverseComplementEdge(e);
                        if (re != null) {
                            reEdgesToRemove.add(re);
                        }                        
                    }
                    graph.removeAllEdges(reEdgesToRemove);
                }
                
                // remove other edges from graph
                graph.removeAllEdges(edgesToRemove);
            }
        }
    }
    
    private ArrayDeque<String> getUnambiguousRightExtension(String id) {
        ArrayDeque<String> path = new ArrayDeque<>();
        HashSet<String> visited = new HashSet<>();
        visited.add(id);
        
        String n = id;
        while (true) {
            List<String> successors = Graphs.successorListOf(graph, n);
            
            if (successors.size() != 1) {
                break;
            }
            else {
                n = successors.get(0);
                List<String> predecessors = Graphs.predecessorListOf(graph, n);
                if (predecessors.size() > 1 || visited.contains(n)) {
                    break;
                }
                else {
                    path.add(n);
                    visited.add(n);
                }
            }
        }
        
        return path;
    }
    
    private ArrayDeque<String> getUnambiguousLeftExtension(String id) {
        ArrayDeque<String> path = new ArrayDeque<>();
        HashSet<String> visited = new HashSet<>();
        visited.add(id);
        
        String n = id;
        while (true) {
            List<String> predecessor = Graphs.predecessorListOf(graph, n);
            
            if (predecessor.size() != 1) {
                break;
            }
            else {
                n = predecessor.get(0);
                List<String> successors = Graphs.successorListOf(graph, n);
                if (successors.size() > 1 || visited.contains(n)) {
                    break;
                }
                else {
                    path.addFirst(n);
                    visited.add(n);
                }
            }
        }
        
        return path;
    }
    
    public void writeBackboneSequences(String outFastaPath) throws IOException {
        if (stranded) {
            layoutStrandedBackbones(outFastaPath);
        }
        else {
            layoutBackbones(outFastaPath);
        }
    }
    
    private String assemblePath(ArrayDeque<String> path, HashMap<String, byte[]> sequences) {
        StringBuilder sb = new StringBuilder();
        
        Iterator<String> itr = path.iterator();
        String vid = itr.next();
        
        boolean reverseComplement = isVertexSignReverseComplement(vid);
        byte[] bytes = sequences.get(getVertexName(vid));
        int start = 0;
        int end = bytes.length;
        
        while (itr.hasNext()) {
            String vid2 = itr.next();
            
            OverlapEdge edge = graph.getEdge(vid, vid2);
            if (reverseComplement) {
                sb.append(reverseComplement(bytes, edge.sourceEnd, end));
            }
            else {
                sb.append(bytesToString(bytes, start, edge.sourceStart));
            }
            
            vid = vid2;
            bytes = sequences.get(getVertexName(vid));
            reverseComplement = isVertexSignReverseComplement(vid);
            
            if (reverseComplement) {
                end = edge.sinkEnd;
                start = 0;
            }
            else {
                start = edge.sinkStart;
                end = bytes.length;
            }
        }
        
        if (reverseComplement) {
            sb.append(reverseComplement(bytes, start, end));
        }
        else {
            sb.append(bytesToString(bytes, start, end));
        }
        
        return sb.toString();
    }
    
    private static class Interval {
        int start;
        int end;
        
        public Interval(int start, int end) {
            this.start = start;
            this.end = end;
        }
    }
    
    private static int getMinCoverage(ArrayDeque<Interval> spans, int targetLength, int edgeTolerance, int window) {
        int effectiveLength = targetLength - 2*edgeTolerance;
        int numWindows = effectiveLength / window;
        int shift = edgeTolerance + (effectiveLength % window)/2;
        
        if (targetLength < window) {
            numWindows = 1;
            shift = 0;
        }
        else if (effectiveLength < window) {
            numWindows = 1;
            shift = (targetLength % window)/2;
        }
        
        int[] covs = new int[numWindows];
        
        for (Interval s : spans) {
            int adjStartPos =  s.start - shift;
            int wStart = adjStartPos/window;
            if (adjStartPos % window > 0) {
                ++wStart;
            }
            
            int wEnd = Math.min((s.end - shift)/window, numWindows);
            
            for (int i=wStart; i<wEnd; ++i) {
                covs[i] += 1;
            }
        }
        
        int min = covs[0];
        for (int c : covs) {
            if (c < min) {
                min = c;
            }
        }
        
        return min;
    }
    
    private void layoutBackbones(String outFastaPath) throws IOException {
        HashMap<String, Integer> lengths = new HashMap<>(); // read id -> length
        HashMap<String, String> longestAlts = new HashMap<>(); // read id -> longest read id
        ArrayDeque<PafRecord> dovetailRecords = new ArrayDeque<>();
        HashMap<String, Integer> artifactCutIndexes = new HashMap<>(); // read id -> cut index
               
        // look for containment and dovetails
        PafReader reader = new PafReader(overlapPafPath);
        
        boolean checkNumAltReads = minNumAltReads > 0;
        ArrayDeque<Interval> spans = new ArrayDeque<>();
        String prevName = null;
        int prevLen = -1;
        ArrayDeque<String> discardReadIDs = new ArrayDeque<>();
        
        while (reader.hasNext()) {
            ExtendedPafRecord r = reader.next();
            
            if (checkNumAltReads) {
                if (!r.qName.equals(prevName)) {
                    if (!spans.isEmpty() && prevName != null) {
                        if (minNumAltReads > getMinCoverage(spans, prevLen, maxEdgeClip, minOverlapMatches)) {
                            discardReadIDs.add(prevName);
                        }
                    }

                    spans = new ArrayDeque<>();
                }

                prevName = r.qName; 
                prevLen = r.qLen;
                spans.add(new Interval(r.qStart, r.qEnd));
            }
            
            lengths.put(r.qName, r.qLen);
            lengths.put(r.tName, r.tLen);
            
            if (r.qName.equals(r.tName)) {
                if (cutRevCompArtifact && 
                        hasReverseComplementArtifact(r) && 
                        hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
                    int cutIndex = getReverseComplementArtifactCutIndex(r);
                    artifactCutIndexes.put(r.qName, cutIndex);
                }
            }
            else {
                if (hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
                    if (isContainmentPafRecord(r)) {
                        String shorter, longer;
                        int longerLen;

                        if (r.qLen > r.tLen || (r.qLen == r.tLen && r.qName.compareTo(r.tName) > 0)) {
                            shorter = r.tName;
                            longer = r.qName;
                            longerLen = r.qLen;
                        }
                        else {
                            shorter = r.qName;
                            longer = r.tName;
                            longerLen = r.tLen;
                        }

                        if (longestAlts.containsKey(shorter)) {
                            String alt = longestAlts.get(shorter);
                            int altLen = lengths.get(alt);
                            if (altLen < longerLen || (altLen == longerLen && longer.compareTo(alt) > 0)) {
                                longestAlts.put(shorter, longer);
                            }
                        }
                        else {
                            longestAlts.put(shorter, longer);
                        }
                    }
                    else if (!longestAlts.containsKey(r.qName) &&
                            !longestAlts.containsKey(r.tName) &&
                            isDovetailPafRecord(r)) {
                        dovetailRecords.add(r);
                    }
                }
            }
        }
        reader.close();
        
        if (checkNumAltReads && !spans.isEmpty() && prevName != null) {
            if (minNumAltReads > getMinCoverage(spans, prevLen, maxEdgeClip, minOverlapMatches)) {
                discardReadIDs.add(prevName);
            }
        }
        
        System.out.println("Overlapped sequences: " + NumberFormat.getInstance().format(lengths.size()));
        if (!discardReadIDs.isEmpty()) {
            System.out.println("         - discarded: " + NumberFormat.getInstance().format(discardReadIDs.size()));
        }
        
        if (!artifactCutIndexes.isEmpty()) {
            System.out.println("         - artifacts: " + NumberFormat.getInstance().format(artifactCutIndexes.size()));
        }
                
        // look for longest reads
        HashSet<String> longestSet = new HashSet<>(longestAlts.values());
        for (String name : lengths.keySet()) {
            if (!longestAlts.containsKey(name) && !longestSet.contains(name)) {
                longestSet.add(name);
            }
        }
        longestSet.removeAll(discardReadIDs);
        
        if (!longestSet.isEmpty()) {
            System.out.println("         - unique:    " + NumberFormat.getInstance().format(longestSet.size()));
        }
        
        // construct overlap graph
        HashSet<String> dovetailReadNames = new HashSet<>(Math.min(longestSet.size(), 2*dovetailRecords.size()));
        for (PafRecord r : dovetailRecords) {
            if ((!cutRevCompArtifact || (!artifactCutIndexes.containsKey(r.qName) && !artifactCutIndexes.containsKey(r.tName))) &&
                    longestSet.contains(r.qName) && longestSet.contains(r.tName)) {
                if (!dovetailReadNames.contains(r.qName)) {
                    graph.addVertex(r.qName + "+");
                    graph.addVertex(r.qName + "-");
                    dovetailReadNames.add(r.qName);
                }
                
                if (!dovetailReadNames.contains(r.tName)) {
                    graph.addVertex(r.tName + "+");
                    graph.addVertex(r.tName + "-");
                    dovetailReadNames.add(r.tName);
                }
                
                if (r.reverseComplemented) {
                    if (r.qEnd >= r.qLen - maxEdgeClip && r.tEnd >= r.tLen - maxEdgeClip) {
                        graph.addEdge(r.tName+"+", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                        graph.addEdge(r.qName+"+", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                    }
                    else if (r.tStart <= maxEdgeClip && r.qStart <= maxEdgeClip) {
                        graph.addEdge(r.qName+"-", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                        graph.addEdge(r.tName+"-", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                    }
                }
                else {
                    if (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) {
                        graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                        graph.addEdge(r.tName+"-", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                    }
                    else if (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip) {
                        graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                        graph.addEdge(r.qName+"-", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                    }
                }
            }
        }
        
        if (!dovetailReadNames.isEmpty()) {
            System.out.println("         - dovetail:  " + NumberFormat.getInstance().format(dovetailReadNames.size()));
        }
        
        int numEdges = graph.edgeSet().size();
//        if (numEdges > 2) {
//            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
//            
//            reduceTransitively();
//            numEdges = graph.edgeSet().size();
//            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
//        }
        
        if (numEdges > 1) {
            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
            //resolveJunctions(dovetailReadNames, false);
            resolveJunctions();
            numEdges = graph.edgeSet().size();
            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
        }
        
        // extract longest read sequences
        HashMap<String, byte[]> longestReadSeqs = new HashMap<>(longestSet.size());
        FastaReader fr = new FastaReader(seqFastaPath);
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        long seqID = 0;
        long originalNumSeq = 0;
        while (fr.hasNext()) {
            ++originalNumSeq;
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            if (dovetailReadNames.contains(name)) {
                longestReadSeqs.put(name, stringToBytes(nameSeq[1], lengths.get(name)));
            }
            else if (longestSet.contains(name) || (!checkNumAltReads && !lengths.containsKey(name))) {
                // an orphan sequence with no overlaps with other sequences
                if (cutRevCompArtifact && artifactCutIndexes.containsKey(name)) {
                    int halfLen = nameSeq[1].length()/2;
                    int cutIndex = artifactCutIndexes.get(name);
                    if (cutIndex < halfLen) {
                        fw.write(Long.toString(++seqID), nameSeq[1].substring(cutIndex+1));
                    }
                    else {
                        fw.write(Long.toString(++seqID), nameSeq[1].substring(0, cutIndex));
                    }
                }
                else {
                    fw.write(Long.toString(++seqID), nameSeq[1]);
                }
            }
        }
        fr.close();
        
        // layout unambiguous paths
        HashSet<String> visitedReadNames = new HashSet<>(dovetailReadNames.size());
        for (String n : dovetailReadNames) {
            if (!visitedReadNames.contains(n)) {
                n += "+";
                
                ArrayDeque<String> path = getUnambiguousLeftExtension(n);
                path.add(n);
                if (!graph.containsEdge(n, path.getFirst())) {
                    // detect cycles
                    path.addAll(getUnambiguousRightExtension(n));
                }
                
                String backbone = assemblePath(path, longestReadSeqs);
                
                // write backbone to FASTA
                String header = Long.toString(++seqID);
                if (path.size() > 1) {
                    //print(path, ',')
                    header += " path=[" + String.join(",", path) + "]";
                }
                fw.write(header, backbone);
                
                for (String vid : path) {
                    visitedReadNames.add(getVertexName(vid));
                }
            }
        }
        fw.close();
        
        System.out.println("before: " + NumberFormat.getInstance().format(originalNumSeq) + "\tafter: " + NumberFormat.getInstance().format(seqID));
    }
    
    private void layoutStrandedBackbones(String outFastaPath) throws IOException {
        HashMap<String, Integer> lengths = new HashMap<>(); // read id -> length
        HashMap<String, String> longestAlts = new HashMap<>(); // read id -> longest read id
        ArrayDeque<PafRecord> dovetailRecords = new ArrayDeque<>();
        HashMap<String, Integer> artifactCutIndexes = new HashMap<>(); // read id -> cut index
        
        // look for containment and overlaps
        PafReader reader = new PafReader(overlapPafPath);
        
        boolean checkNumAltReads = minNumAltReads > 0;
        ArrayDeque<Interval> spans = new ArrayDeque<>();
        String prevName = null;
        int prevLen = -1;
        ArrayDeque<String> discardReadIDs = new ArrayDeque<>();
        
        while (reader.hasNext()) {
            ExtendedPafRecord r = reader.next();

            if (checkNumAltReads && !r.reverseComplemented) {
                if (!r.qName.equals(prevName)) {
                    if (!spans.isEmpty() && prevName != null) {
                        if (minNumAltReads > getMinCoverage(spans, prevLen, maxEdgeClip, minOverlapMatches)) {
                            discardReadIDs.add(prevName);
                        }
                    }

                    spans = new ArrayDeque<>();
                }

                prevName = r.qName; 
                prevLen = r.qLen;
                spans.add(new Interval(r.qStart, r.qEnd));
            }
            
            lengths.put(r.qName, r.qLen);
            lengths.put(r.tName, r.tLen);
            
            if (r.qName.equals(r.tName)) {
                if (cutRevCompArtifact && 
                        hasReverseComplementArtifact(r) && 
                        hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
                    int cutIndex = getReverseComplementArtifactCutIndex(r);
                    artifactCutIndexes.put(r.qName, cutIndex);
                }
            }
            else {
                if (hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
                    if (isStrandedContainmentPafRecord(r)) {
                        String shorter, longer;
                        int longerLen;

                        if (r.qLen > r.tLen || (r.qLen == r.tLen && r.qName.compareTo(r.tName) > 0)) {
                            shorter = r.tName;
                            longer = r.qName;
                            longerLen = r.qLen;
                        } 
                        else {
                            shorter = r.qName;
                            longer = r.tName;
                            longerLen = r.tLen;
                        }

                        if (longestAlts.containsKey(shorter)) {
                            String alt = longestAlts.get(shorter);
                            int altLen = lengths.get(alt);
                            if (altLen < longerLen || (altLen == longerLen && longer.compareTo(alt) > 0)) {
                                longestAlts.put(shorter, longer);
                            }
                        }
                        else {
                            longestAlts.put(shorter, longer);
                        }
                    }
                    else if (!longestAlts.containsKey(r.qName) &&
                            !longestAlts.containsKey(r.tName) &&
                            isStrandedDovetailPafRecord(r)) {
                        dovetailRecords.add(r);
                    }
                }
            }
        }
        reader.close();
        
        if (checkNumAltReads && !spans.isEmpty() && prevName != null) {
            if (minNumAltReads > getMinCoverage(spans, prevLen, maxEdgeClip, minOverlapMatches)) {
                discardReadIDs.add(prevName);
            }
        }
        
        System.out.println("Overlapped sequences: " + NumberFormat.getInstance().format(lengths.size()));
        if (!discardReadIDs.isEmpty()) {
            System.out.println("         - discarded: " + NumberFormat.getInstance().format(discardReadIDs.size()));
        }
        
        if (!artifactCutIndexes.isEmpty()) {
            System.out.println("         - artifacts: " + NumberFormat.getInstance().format(artifactCutIndexes.size()));
        }
        
        // look for longest reads
        HashSet<String> longestSet = new HashSet<>(longestAlts.values());
        for (String name : lengths.keySet()) {
            if (!longestAlts.containsKey(name) && !longestSet.contains(name)) {
                longestSet.add(name);
            }
        }
        longestSet.removeAll(discardReadIDs);
        
        if (!longestSet.isEmpty()) {
            System.out.println("         - unique:    " + NumberFormat.getInstance().format(longestSet.size()));
        }        
        
        // construct overlap graph
        HashSet<String> dovetailReadNames = new HashSet<>(Math.min(longestSet.size(), 2*dovetailRecords.size()));
        for (PafRecord r : dovetailRecords) {
            if (!r.reverseComplemented) {
                if ((!cutRevCompArtifact || (!artifactCutIndexes.containsKey(r.qName) && !artifactCutIndexes.containsKey(r.tName))) &&
                        longestSet.contains(r.qName) && longestSet.contains(r.tName)) {
                    if (!dovetailReadNames.contains(r.qName)) {
                        graph.addVertex(r.qName + "+");
                        dovetailReadNames.add(r.qName);
                    }
                    
                    if (!dovetailReadNames.contains(r.tName)) {
                        graph.addVertex(r.tName + "+");
                        dovetailReadNames.add(r.tName);
                    }
                    
                    if (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) {
                        graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                    }
                    else if (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip) {
                        graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                    }
                }
            }
        }
        
        if (!dovetailReadNames.isEmpty()) {
            System.out.println("         - dovetail:  " + NumberFormat.getInstance().format(dovetailReadNames.size()));
        }
        
        int numEdges = graph.edgeSet().size();
//        if (numEdges > 2) {
//            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
//            
//            reduceTransitively();
//            numEdges = graph.edgeSet().size();
//            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
//        }
        
        if (numEdges > 1) {
            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
            //resolveJunctions(dovetailReadNames, true);
            resolveJunctions();
            numEdges = graph.edgeSet().size();
            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
        }
        
        // extract longest read sequences
        HashMap<String, byte[]> longestReadSeqs = new HashMap<>(dovetailReadNames.size());
        FastaReader fr = new FastaReader(seqFastaPath);
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        long seqID = 0;
        long originalNumSeq = 0;
        while (fr.hasNext()) {
            ++originalNumSeq;
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            if (dovetailReadNames.contains(name)) {
                String seq = nameSeq[1];
                longestReadSeqs.put(name, stringToBytes(seq, seq.length()));
            }
            else if (longestSet.contains(name) || (!checkNumAltReads && !lengths.containsKey(name))) {
                // an orphan sequence with no overlaps with other sequences
                if (cutRevCompArtifact && artifactCutIndexes.containsKey(name)) {
                    int halfLen = nameSeq[1].length()/2;
                    int cutIndex = artifactCutIndexes.get(name);
                    if (cutIndex < halfLen) {
                        fw.write(Long.toString(++seqID), nameSeq[1].substring(cutIndex+1));
                    }
                    else {
                        fw.write(Long.toString(++seqID), nameSeq[1].substring(0, cutIndex));
                    }
                }
                else {
                    fw.write(Long.toString(++seqID), nameSeq[1]);
                }
            }
        }
        fr.close();
        
        // layout unambiguous paths
        HashSet<String> visitedReadNames = new HashSet<>(dovetailReadNames.size());
        for (String n : dovetailReadNames) {
            if (!visitedReadNames.contains(n)) {
                n += "+";
                
                ArrayDeque<String> path = getUnambiguousLeftExtension(n);
                path.add(n);
                if (!graph.containsEdge(n, path.getFirst())) {
                    // detect cycles
                    path.addAll(getUnambiguousRightExtension(n));
                }
                
                //print(path, ',')
                String backbone = assemblePath(path, longestReadSeqs);
                
                // write backbone to FASTA
                String header = Long.toString(++seqID);
                if (path.size() > 1) {
                    //print(path, ',')
                    header += " path=[" + String.join(",", path) + "]";
                }
                fw.write(header, backbone);
                
                for (String vid : path) {
                    visitedReadNames.add(getVertexName(vid));
                }
            }
        }
        fw.close();
        
        System.out.println("before: " + NumberFormat.getInstance().format(originalNumSeq) + "\tafter: " + NumberFormat.getInstance().format(seqID));
    }
    
    public static void main(String[] args) {
        boolean stranded = true;
        String seqFastaPath = "cat.fa";
        String overlapPafPath = "ava.paf.gz";
        String backboneFastaPath = "backbones.fa";
        
        try {
            //seqFile, pafFile, stranded, maxEdgeClip, minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact
            Layout myLayout = new Layout(seqFastaPath, overlapPafPath, stranded, 5, 0.90f, 2*25, 1, true, 1);
            myLayout.writeBackboneSequences(backboneFastaPath);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
