/* 
 * Copyright (C) 2018-present BC Cancer Genome Sciences Centre
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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.Writer;
import java.text.NumberFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import static org.jgrapht.Graph.DEFAULT_EDGE_WEIGHT;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultEdge;
import rnabloom.util.MiniFloat;
import rnabloom.io.CompressedFastaRecord;
import static rnabloom.io.Constants.FASTA_EXT;
import static rnabloom.io.Constants.GZIP_EXT;
import rnabloom.io.ExtendedPafRecord;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaRecord;
import rnabloom.io.FastaWriter;
import rnabloom.io.PafReader;
import rnabloom.io.PafRecord;
import rnabloom.util.BitSequence;
import static rnabloom.util.Common.convertToRoundedPercent;
import static rnabloom.util.FileUtils.fileToStringCollection;
import static rnabloom.util.IntervalUtils.getOverlap;
import static rnabloom.util.IntervalUtils.isContained;
import static rnabloom.util.IntervalUtils.merge;
import static rnabloom.util.PafUtils.*;
import rnabloom.util.PolyATailFinder;
import static rnabloom.util.SeqUtils.isLowComplexityLongWindowed;
import rnabloom.util.Timer;
import static rnabloom.util.FileUtils.getTextFileWriter;
import static rnabloom.util.FileUtils.readIntArrayFromFile;
import static rnabloom.util.IntervalUtils.isForwardDoveTail;
import smile.stat.distribution.BinomialDistribution;
import smile.stat.distribution.EmpiricalDistribution;
import smile.util.IntSet;

/**
 *
 * @author Ka Ming Nip
 */
public class Layout {
    
    private DefaultDirectedWeightedGraph<String, OverlapEdge> graph;
    private InputStream overlapPafInputStream;
    private String seqFastaPath;
    private boolean stranded;
    private int maxEdgeClip = 100;
    private float minAlnId = 0.50f;
    private int maxIndelSize = 20;
    private int minOverlapMatches = 200;
    private boolean cutRevCompArtifact = false;
    private int minNumAltReads = 0;
    private boolean verbose = false;
    HashMap<String, PolyAInfo> polyAInfoMap = null;
    private int minPolyALength = 0;
    private final OverlapSizeComparator overlapSizeComparator = new OverlapSizeComparator();
    private final OverlapQueryStartComparator qStartComparator = new OverlapQueryStartComparator();
    
    public Layout(String seqFile, InputStream overlapPafInputStream, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth, boolean verbose,
            PolyATailFinder polyaFinder) throws IOException {
        this.graph = new DefaultDirectedWeightedGraph<>(OverlapEdge.class);
        this.overlapPafInputStream = overlapPafInputStream;
        this.seqFastaPath = seqFile;
        this.stranded = stranded;
        this.maxEdgeClip = maxEdgeClip;
        this.minAlnId = minAlnId;
        this.minOverlapMatches = minOverlapMatches;
        this.maxIndelSize = maxIndelSize;
        this.cutRevCompArtifact = cutRevCompArtifact;
        this.minNumAltReads = minSeqDepth - 1;
        this.verbose = verbose;
        if (polyaFinder != null) {
            this.minPolyALength = polyaFinder.getSeedLength();
            this.polyAInfoMap = getPolyAInfo(polyaFinder);
        }
    }
    
    private void printMessage(String msg) {
        if (verbose) {
            System.out.println(msg);
        }
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
            //return Math.max(o.sinkEnd-o.sinkStart, o.sourceEnd-o.sourceStart) - Math.max(sinkEnd-sinkStart, sourceEnd-sourceStart);
            return ((o.sinkEnd-o.sinkStart)+(o.sourceEnd-o.sourceStart))/2 - ((sinkEnd-sinkStart)+(sourceEnd-sourceStart))/2;
        }
    }

    private boolean hasLargeOverlap(PafRecord r) {
        return rnabloom.util.PafUtils.hasLargeOverlap(r, minOverlapMatches);
    }
    
    private boolean hasGoodOverlap(PafRecord r) {
        return rnabloom.util.PafUtils.hasGoodOverlap(r, minAlnId);
    }
        
    private boolean hasGoodAlignment(ExtendedPafRecord r) {
        return rnabloom.util.PafUtils.hasGoodAlignment(r, maxIndelSize, minAlnId);
    }
    
    private boolean hasReverseComplementArtifact(ExtendedPafRecord r) {
        return rnabloom.util.PafUtils.hasReverseComplementArtifact(r, maxEdgeClip);
    }
    
    private int getReverseComplementArtifactCutIndex(ExtendedPafRecord r) {
        if (r.qName.equals(r.tName) && r.reverseComplemented) {
            int rightEdgeClip = r.qLen - r.qEnd;
            
            if (r.qStart <= maxEdgeClip && r.qStart <= rightEdgeClip) {
                // left edge has reverse complement
                if (r.qStart == r.tStart && r.qEnd == r.tEnd) {
                    // palindrome
                    return r.qStart + (r.qEnd - r.qStart)/2;
                }
                else {
                    return r.qEnd -1;
                }
            }
            else if (rightEdgeClip <= maxEdgeClip && r.qStart >= rightEdgeClip) {
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
    
    private boolean isContainmentPafRecord(PafRecord r) {
        return rnabloom.util.PafUtils.isContainmentPafRecord(r, maxEdgeClip);
    }
    
    private boolean isForwardContainmentPafRecord(PafRecord r) {
        return !r.reverseComplemented && isContainmentPafRecord(r);
    }
        
    private boolean isDovetailPafRecord(PafRecord r) {
        return rnabloom.util.PafUtils.isDovetailPafRecord(r, maxEdgeClip);
    }
    
    private boolean isForwardDovetailPafRecord(PafRecord r) { 
        return rnabloom.util.PafUtils.isForwardDovetailPafRecord(r, maxEdgeClip);
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
    
    private void removeTransitiveEdges() {
        HashSet<String> walkedSet = new HashSet<>(graph.vertexSet().size());
        
        ArrayDeque<OverlapEdge> edgesToRemove = new ArrayDeque<>();
        
        for (String vid : graph.vertexSet()) {
            if (!walkedSet.contains(vid)) {
                // greedy acyclic walk
                ArrayDeque<String> path = getGreedyExtension(vid);

                // find transitive edges on the path
                HashSet<String> pathSet = new HashSet<>(path);
                String prev = null;
                for (String cursor : path) {
                    if (prev != null) {
                        Set<OverlapEdge> outEdges = graph.outgoingEdgesOf(prev);
                        if (outEdges.size() > 1) {
                            for (OverlapEdge e : outEdges) {
                                String target = graph.getEdgeTarget(e);
                                if (pathSet.contains(target) && !target.equals(cursor)) {
                                    edgesToRemove.add(e);
                                }
                            }
                        }
                    }

                    prev = cursor;
                }

                // remove transitive edges
                graph.removeAllEdges(edgesToRemove);
                edgesToRemove.clear();
                
                // store path vertexes
                walkedSet.addAll(path);
            }
        }
    }
    
    private ArrayDeque removeRedundantNodes() {
        ArrayDeque<String> removed = new ArrayDeque<>();
        
        for (String name : new ArrayDeque<>(graph.vertexSet())) {
            if (isRedundantNode(name)) {
                removeVertexStranded(name);
                removed.add(name);
            }
        }
        
        return removed;
    }
    
    private boolean isRedundantNode(String name) {
        Set<OverlapEdge> inEdgeSet = graph.incomingEdgesOf(name);
        if (inEdgeSet.isEmpty()) {
            // this node is a leaf
            return false;
        }
        
        Set<OverlapEdge> outEdgeSet = graph.outgoingEdgesOf(name);
        if (outEdgeSet.isEmpty()) {
            // this node is a leaf
            return false;
        }
        
        ArrayList<OverlapEdge> inEdges = new ArrayList<>(inEdgeSet);
        ArrayList<OverlapEdge> outEdges = new ArrayList<>(outEdgeSet);
        
        Collections.sort(inEdges);
        Collections.sort(outEdges);
                
        String closetPredecessor = graph.getEdgeSource(inEdges.get(0));
        String closetSuccessor = graph.getEdgeTarget(outEdges.get(0));
        
        OverlapEdge bridgeEdge = graph.getEdge(closetPredecessor, closetSuccessor);
        
        if (bridgeEdge == null) {
            // this node is not redundant
            return false;
        }
        
        int numPredecessors = inEdges.size();
        int numSuccessors = outEdges.size();
        
        HashSet<String> predecessorSet = new HashSet<>(numPredecessors);
        for (OverlapEdge e : inEdges) {
            String p = graph.getEdgeSource(e);
            predecessorSet.add(p);
        }
        
        HashSet<String> successorSet = new HashSet<>(numSuccessors);
        for (OverlapEdge e : outEdges) {
            String s = graph.getEdgeTarget(e);
            successorSet.add(s);
        }
                
        // each predecessor must have one outgoing edges leading to another predecessor/successor
        HashSet<String> pendingPredecessors = new HashSet<>();
        HashSet<String> bridgedPredecessors = new HashSet<>();
        HashSet<String> bridgedSuccessors = new HashSet<>();
        for (String p : predecessorSet) {
            Set<OverlapEdge> edges = graph.outgoingEdgesOf(p);
            OverlapEdge inEdge = graph.getEdge(p, name);
            boolean foundBridge = false;
            for (OverlapEdge e : edges) {
                if (e != inEdge) {
                    String s = graph.getEdgeTarget(e);
                    if (successorSet.contains(s)) {
                        // make sure distance is similar
                        int d = ((e.sinkEnd - e.sinkStart) + (e.sourceEnd - e.sourceStart))/2;
                        
                        OverlapEdge outEdge = graph.getEdge(name, s);
                        int length = isVertexSignReverseComplement(name) ?
                                inEdge.sinkEnd - outEdge.sourceStart : 
                                outEdge.sourceEnd - inEdge.sinkStart;                        
                        int inEdgeNotCovered = length - (inEdge.sinkEnd - inEdge.sinkStart);
                        int outEdgeNotCovered = length - (outEdge.sourceEnd - outEdge.sourceStart);
                        int d2 = length - inEdgeNotCovered - outEdgeNotCovered;
                        
                        if (Math.max(d, d2) * 0.9f > Math.min(d, d2)) {
                            return false;
                        }
                        else {
                            foundBridge = true;
                            bridgedSuccessors.add(s);
                        }
                    }
                }
            }
            
            if (foundBridge) {
                bridgedPredecessors.add(p);
            }
            else {
                pendingPredecessors.add(p);
            }
        }
        
        // check pending predecessors
        if (!pendingPredecessors.isEmpty()) {
            for (String p : pendingPredecessors) {
                boolean found = false;
                for (String s : Graphs.successorListOf(graph, p)) {
                    if (bridgedPredecessors.contains(s)) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    return false;
                }
            }
        }
        
        // check pending successors
        successorSet.removeAll(bridgedSuccessors);
        if (!successorSet.isEmpty()) {
            for (String s : successorSet) {
                boolean found = false;
                for (String p : Graphs.predecessorListOf(graph, s)) {
                    if (bridgedSuccessors.contains(p)) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    return false;
                }
            }
        }
        
        return true;
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
    
    private ArrayDeque<String> getUnambiguousExtension(String vid) {
        ArrayDeque<String> path = new ArrayDeque<>();
        HashSet<String> visited = new HashSet<>();
        visited.add(vid);
        
        // extend left
        String n = vid;
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
        
        path.add(vid);
        
        // extend right
        n = vid;
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
        
    private String assemblePath(ArrayDeque<String> path, HashMap<String, BitSequence> sequences) {
        StringBuilder sb = new StringBuilder();
        
        Iterator<String> itr = path.iterator();
        String vid = itr.next();
        
        boolean reverseComplement = isVertexSignReverseComplement(vid);
        BitSequence bitseq = sequences.get(getVertexName(vid));
        int start = 0;
        int end = bitseq.length;
        
        while (itr.hasNext()) {
            String vid2 = itr.next();
            
            OverlapEdge edge = graph.getEdge(vid, vid2);            
            if (reverseComplement) {
                if (end > edge.sourceEnd) {
                    sb.append(bitseq.subStringRevComp(edge.sourceEnd, end));
                }
            }
            else {
                if (edge.sourceStart > start) {
                    sb.append(bitseq.subString(start, edge.sourceStart));
                }
            }
            
            vid = vid2;
            bitseq = sequences.get(getVertexName(vid));
            reverseComplement = isVertexSignReverseComplement(vid);
            
            if (reverseComplement) {
                end = edge.sinkEnd;
                start = 0;
            }
            else {
                start = edge.sinkStart;
                end = bitseq.length;
            }
        }
        
        if (reverseComplement) {
            sb.append(bitseq.subStringRevComp(start, end));
        }
        else {
            sb.append(bitseq.subString(start, end));
        }
        
        return sb.toString();
    }

    /*
    private static Interval overlapInterval(Interval span1, Interval span2) {
        if (span1.start >= span2.start && span1.start <= span2.end) {
            return new Interval(span2.start, Math.max(span1.end, span2.end));
        }
        else if (span1.end >= span2.start && span1.end <= span2.end) {
            return new Interval(Math.min(span1.start, span2.start), span2.end);
        }
                
        return null;
    }
    
    private static ArrayList<Interval> overlapIntervals(ArrayList<Interval> spans) {
        int numSpans = spans.size();
        
        for (int cursor=0; cursor<numSpans; ++cursor) {
            ArrayList<Interval> newSpans = new ArrayList<>();
            int newCursor = 0;
            Interval query = spans.get(cursor);

            for (int i=0; i<cursor; ++i) {
                Interval target = spans.get(i);
                Interval merged = overlapInterval(query, target);
                if (merged != null) {
                    query = merged;
                }
                else {
                    newSpans.add(target);
                }
            }
            
            newCursor = newSpans.size();
            newSpans.add(query);

            for (int i=cursor+1; i<numSpans; ++i) {
                Interval target = spans.get(i);
                Interval merged = overlapInterval(query, target);
                if (merged != null) {
                    query = merged;
                    newSpans.set(newCursor, merged);
                }
                else {
                    newSpans.add(target);
                }
            }
            
            spans = newSpans;
            numSpans = newSpans.size();
            cursor = newCursor;
        }
        
        return spans;
    }
    */
    
    private static void overlapIntervals(TreeSet<ComparableInterval> targets) {
        if (!targets.isEmpty()) {
            Iterator<ComparableInterval> itr = targets.iterator();
            ComparableInterval prev = itr.next();
            while (itr.hasNext()) {
                ComparableInterval i = itr.next();
                if (prev.merge(i)) {
                    itr.remove();
                }
                else {
                    prev = i;
                }
            }
        }
    }
    
    private static void overlapIntervals(TreeSet<ComparableInterval> targets, ComparableInterval query) {
        Iterator<ComparableInterval> itr = targets.iterator();
        while (itr.hasNext()) {
            ComparableInterval i = itr.next();
            if (query.merge(i)) {
                itr.remove();
            }
            else if (query.compareTo(i) < 0) {
                break;
            }
        }
        targets.add(query);
    }
        
    private static ArrayDeque<Interval> extractEffectiveIntervals(Collection<Interval> spans, int minCoverage, int minIntervalLength) {
        ArrayDeque<Interval> effIntervals = new ArrayDeque<>();
        
        int targetLength = 0;
        for (Interval s : spans) {
            if (s.end > targetLength) {
                targetLength = s.end;
            }
        }
        
        // extract coverage histogram
        int[] hist = new int[targetLength];
        for (Interval s : spans) {
            for (int i=s.start; i<s.end; ++i) {
                ++hist[i];
            }
        }
        
        // extract effective intervals
        int start = -1;
        for (int i=0; i<targetLength; ++i) {
            int c = hist[i];
            if (c >= minCoverage) {
                if (start < 0) {
                    start = i;
                }
            }
            else {
                if (start >= 0) {
                    if (minIntervalLength <= i - start) {
                        effIntervals.add(new Interval(start, i));
                    }
                    start = -1;
                }
            }
        }

        // extract last interval
        if (start >= 0 && minIntervalLength <= targetLength - start) {
            effIntervals.add(new Interval(start, targetLength));
        }
        
        return effIntervals;
    }
    
    private static int getMinCoverage(Collection<Interval> spans, int targetLength, int edgeTolerance, int window) {
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
    
    private int getHistogramBinSize(int length) {
        /* Set binsize to `minOverlapMatches/2` to ensure that 2 consecutive bins
            represent reads overlapping by >= `minOverlapMatches` because only
            the middle bin(s) of each interval will contribute to the histogram
                    |---|===|---|
                |---|===|---|
            |---|===|---|
                 === === ===
            | 0 | 1 | 1 | 1 | 0 |
        */
        if (length <= 250) {
            return Math.min(25, minOverlapMatches/2);
        }
        else if (length <= 500) {
            return Math.min(50, minOverlapMatches/2);
        }
        else if (length <= 1000) {
            return Math.min(100, minOverlapMatches/2);
        }
        else {
            return Math.min(200, minOverlapMatches/2);
        }
    }
    
    private class Histogram {
        int length, minStart, maxEnd = -1;
        byte[] bars = null;
        ArrayDeque<QueryOverlap> pendingQueries = null;
        boolean seenAsQuery = false;
        
        public Histogram(int length, int minStart, int maxEnd, int binSize) {
            this.length = length;
            this.minStart = minStart;
            this.maxEnd = maxEnd;
            this.bars = initByteHistogram(length, binSize);
        }
        
    }
    
    private static byte[] initByteHistogram(int origLength, int binSize) {
        return new byte[(int)Math.ceil((float)origLength/(float)binSize)];
    }

    private static short[] initShortHistogram(int origLength, int binSize) {
        return new short[(int)Math.ceil((float)origLength/(float)binSize)];
    }

    private static int[] initIntHistogram(int origLength, int binSize) {
        return new int[(int)Math.ceil((float)origLength/(float)binSize)];
    }
    
    private static void updateHistogram(Histogram h, int start, int end, int binSize) {
        h.minStart = Math.min(h.minStart, start);
        h.maxEnd = Math.max(h.maxEnd, end);

        if (h.bars != null) {
            updateHistogram(h.bars, h.length, start, end, binSize);
        }
    }
    
    private static void updateHistogram(byte[] h, int length, int start, int end, int binSize) {
        int numBars = h.length;

        if (start > 0) {
            start = (int) Math.round(start/(float)binSize) + 1;
            // +1 to ensure that stacking reads overlap
        }
        else {
            start = 0;
        }

        if (end < length) {
            end = (int) Math.round(end/(float)binSize) - 1;
            // -1 to ensure that stacking reads overlap
        }
        else {
            end = numBars;
        }

        if (start < end && start >= 0 && start < numBars && end > 0 && end <= numBars) {
            updateHistogram(h, start, end);
        }
    }
    
    private boolean isFullyCovered(Histogram h) {
        //return h.minStart <= maxEdgeClip && h.maxEnd >= h.length - maxEdgeClip;
        return h.minStart <= 1 && h.maxEnd >= h.length - 1;
    }
    
    private static int updateHistogram(int[] hist, int binSize, int origLen, int start, int end) {
        int numBars = hist.length;
        
        if (start > 0) {
            start = (int) Math.ceil((float) start/binSize) + 1;
        }
        else {
            start = 0;
        }
        
        if (end < origLen) {
            end = (int) Math.floor((float) end/binSize) - 1;
        }
        else {
            end = numBars;
        }
        
        if (start < end && start >= 0 && start < numBars && end > 0 && end <= numBars) {
            // update histogram and return the min updated count in the interval
            int min = Integer.MAX_VALUE;
            
            for (int i=start; i<end; ++i) {
                hist[i] += 1;
                min = Math.min(min, hist[i]);
            }

            if (min < Integer.MAX_VALUE) {
                return min;
            }
        }
        
        return 0;
    }
    
    private static void updateHistogram(short[] hist, int binSize, int origLen, int start, int end) {
        if (start > 0) {
            start = (int) Math.ceil((float) start/binSize);
        }
        else {
            start = 0;
        }
        
        if (end < origLen) {
            end = (int) Math.floor((float) end/binSize);
        }
        else {
            end = hist.length;
        }
        
        updateHistogram(hist, start, end);
    }
    
    private static void updateHistogram(short[] hist, int start, int end) {
        short c;
        for (int i=start; i<end; ++i) {
            c = hist[i];
            if (c < Short.MAX_VALUE) {
                hist[i] = (short) (c + 1);
            }
        }
    }

    private static void updateHistogram(byte[] hist, int start, int end) {
        for (int i=start; i<end; ++i) {
            hist[i] = MiniFloat.increment(hist[i]);
        }
    }
    
    private static boolean isMultiSegmentHistogram(short[] hist, int minCoverage, int minSegmentLength) {
        int histLength = hist.length;
        int numSegments = 0;
        
        // extract effective intervals
        int start = -1;
        for (int i=0; i<histLength; ++i) {
            int c = hist[i];
            if (c >= minCoverage) {
                if (start < 0) {
                    start = i;
                }
            }
            else {
                if (start >= 0) {
                    if (minSegmentLength <= i - start) {
                        if (++numSegments >= 2) {
                            return true;
                        }
                    }
                    start = -1;
                }
            }
        }

        // extract last interval
        if (start >= 0 && minSegmentLength <= histLength - start) {
            ++numSegments;
        }
    
        return numSegments >= 2;
    }
    
    private static ArrayDeque<Interval> extractEffectiveIntervals(Histogram hist, int binSize, int minCoverage, int minIntervalLength) {
        ArrayDeque<Interval> intervals = extractEffectiveIntervals(hist.bars, binSize, minCoverage, minIntervalLength);
        
        if (!intervals.isEmpty()) {
            Interval first = intervals.getFirst();
            Interval last = intervals.getLast();
            
            /*
            Need to adjust the first and last intervals here because both edges of
            the histogram would not be populated if the alignments do not precisely
            start at 0 or terminate at the end of the sequence.
            */
            int buffer = binSize;
            
            if (first.start < hist.minStart + buffer) {
                first.start = hist.minStart;
            }

            if (last.end > hist.length) {
                last.end = hist.length;
            }
            else if (last.end > hist.maxEnd - buffer) {
                last.end = hist.maxEnd;
            }
        }
        
        return intervals;
    }
    
    private static ArrayDeque<Interval> extractEffectiveIntervals(short[] hist, int binSize, int minCoverage, int minIntervalLength) {
        ArrayDeque<Interval> effIntervals = new ArrayDeque<>();
        int minHistIntervalLength = (int) Math.floor((float)minIntervalLength/(float)binSize);
        int histLength = hist.length;
        
        // extract effective intervals
        int start = -1;
        for (int i=0; i<histLength; ++i) {
            int c = hist[i];
            if (c >= minCoverage) {
                if (start < 0) {
                    start = i;
                }
            }
            else {
                if (start >= 0) {
                    if (minHistIntervalLength <= i - start) {
                        effIntervals.add(new Interval(start*binSize, i*binSize));
                    }
                    start = -1;
                }
            }
        }

        // extract last interval
        if (start >= 0 && minHistIntervalLength <= histLength - start) {
            // histLength*binSize can be > actual sequence length
            effIntervals.add(new Interval(start*binSize, histLength*binSize));
        }
        
        return effIntervals;
    }
    
    private static ArrayDeque<Interval> extractEffectiveIntervals(byte[] hist, int binSize, int minCoverage, int minIntervalLength) {
        ArrayDeque<Interval> effIntervals = new ArrayDeque<>();
        int minHistIntervalLength = (int) Math.floor((float)minIntervalLength/(float)binSize);
        int histLength = hist.length;
        
        // extract effective intervals
        int start = -1;
        for (int i=0; i<histLength; ++i) {
            float c = MiniFloat.toFloat(hist[i]);
            if (c >= minCoverage) {
                if (start < 0) {
                    start = i;
                }
            }
            else {
                if (start >= 0) {
                    int length = i - start;
                    int end = i;
                    
                    /* increase the interval by one on both ends because the first and last
                    histogram bars for each alignment are not incremented to ensure the
                    minimum overlap size between reads.*/
                    if (start > 1) {
                        --start;
                        ++length;
                    }
                    if (i <= histLength-2) {
                        ++end;
                        ++length;
                    }
                    
                    if (minHistIntervalLength <= length) {
                        effIntervals.add(new Interval(start*binSize, end*binSize));
                    }
                    start = -1;
                }
            }
        }

        // extract last interval
        if (start >= 0 && minHistIntervalLength <= histLength - start) {
            // histLength*binSize can be > actual sequence length
            effIntervals.add(new Interval(start*binSize, histLength*binSize));
        }
        
        return effIntervals;
    }
    
    private class ReadGroup {
        public int id = -1;
        private HashSet<String> names = new HashSet<>();
    }
    
    private class ReadClusters3 {
        private int numClusters = 0;
        private HashMap<String, ReadGroup> assignment = new HashMap<>();
        private int mergedClusterMaxSize = 10000;
        private int largestClusterSize = 0;
        private int largestClusterID = -1;
        
        public ReadClusters3() {
            
        }
        
        public ReadClusters3(int mergedClusterMaxSize) {
            this.mergedClusterMaxSize = mergedClusterMaxSize;
        }
        
        public void add(String name1, String name2){
            ReadGroup cluster1 = assignment.get(name1);
            ReadGroup cluster2 = assignment.get(name2);
            
            if (cluster1 == null && cluster2 == null) {
                ReadGroup cluster = new ReadGroup();
                cluster.names.add(name1);
                cluster.names.add(name2);
                assignment.put(name1, cluster);
                assignment.put(name2, cluster);
                ++numClusters;
            }
            else if (cluster1 == null) {
                cluster2.names.add(name1);
                assignment.put(name1, cluster2);
            }
            else if (cluster2 == null) {
                cluster1.names.add(name2);
                assignment.put(name2, cluster1);
            }
            else if (cluster1 != cluster2) {
                int clusterSize1 = cluster1.names.size();
                int clusterSize2 = cluster2.names.size();
                if (clusterSize1 + clusterSize2 < mergedClusterMaxSize) {
                    // merge clusters
                    
                    ReadGroup larger = cluster1;
                    ReadGroup smaller = cluster2;
                    
                    if (clusterSize2 > clusterSize1) {
                        larger = cluster2;
                        smaller = cluster1;
                    }
                    
                    larger.names.addAll(smaller.names);
                    
                    for (String r : smaller.names) {
                        assignment.put(r, larger);
                    }
                    
                    --numClusters;
                }
            }
        }
        
        public void assignIDs() {                        
            int cid = 0;
            for (ReadGroup c : assignment.values()) {
                if (c.id < 0) {
                    c.id = ++cid;
                    
                    int clusterSize = c.names.size();
                    if (clusterSize > largestClusterSize) {
                        largestClusterSize = clusterSize;
                        largestClusterID = cid;
                    }
                }
            }
            
            numClusters = cid;
        }
        
        public int getID(String name) {
            ReadGroup r = assignment.get(name);
            if (r == null) {
                return -1;
            }
            return r.id;
        }
        
        public int size() {
            return numClusters;
        }
        
        public int getLargestClusterSize() {
            return largestClusterSize;
        }
        
        public int getLargestClusterID() {
            return largestClusterID;
        }
    }
    
    private class NeighborPair implements Comparable<NeighborPair> {
        String name1;
        String name2;
        int score;
        
        public NeighborPair(String name1, String name2, int score) {
            this.name1 = name1;
            this.name2 = name2;
            this.score = score;
        }
        
        public boolean equals(NeighborPair o) {
            return this.score == o.score && 
                    name1.equals(o.name1) && name2.equals(o.name2);
        }
                
        @Override
        public int compareTo(NeighborPair o) {
            int diff = o.score - this.score;
            if (diff != 0) {
                return diff;
            }
            else if (!equals(o)) {
                diff = o.name1.compareTo(name1);
                if (diff != 0) {
                    return diff;
                }
                else {
                    return o.name2.compareTo(name2);
                }
            }
            else {
                return 0;
            }
        }
    }
    
    private class Neighbor implements Comparable<Neighbor> {
        String name;
        int score;
        
        public Neighbor(String name, int score) {
            this.name = name;
            this.score = score;
        }
        
        @Override
        public int compareTo(Neighbor o) {
            return o.score - this.score; 
        }
    }
    
    private class BestNeighbor {
        public HashMap<String, Neighbor> neighbors = new HashMap<>();
        
        public BestNeighbor() {
        }
        
        public void add(String name1, String name2, int score) {
            addHelper(name1, name2, score);
            addHelper(name2, name1, score);
        }
        
        private void addHelper(String target, String query, int score) {
            Neighbor n = neighbors.get(target);
            if (n == null) {
                neighbors.put(target, new Neighbor(query, score));
            }
            else if (score > n.score) {
                n.name = query;
                n.score = score;
            }
        }
    }
    
    private class BestNeighborPairs {
        public int max = 2;
        private int maxIndex = max-1; 
        private HashMap<String, ArrayList<NeighborPair>> neighbors = new HashMap<>();
        
        public BestNeighborPairs() {
        }
        
        public BestNeighborPairs(int max) {
            this.max = max;
            this.maxIndex = max - 1;
        }
        
        public void add(String name1, String name2, int score) {
            boolean add1 = false;
            boolean add2 = false;
            
            ArrayList<NeighborPair> arr1 = neighbors.get(name1);
            if (arr1 != null) {
                if (arr1.size() < max) {
                    add1 = true;
                }
                else {
                    NeighborPair worst = arr1.get(maxIndex);
                    if (score > worst.score) {
                        arr1.remove(max-1);
                        add1 = true;
                    }
                }
            }
            else {
                arr1 = new ArrayList<>(max);
                add1 = true;
                neighbors.put(name1, arr1);
            }
            
            ArrayList<NeighborPair> arr2 = neighbors.get(name2);
            if (arr2 != null) {
                if (arr2.size() < max) {
                    add2 = true;
                }
                else {
                    NeighborPair worst = arr2.get(maxIndex);
                    if (score > worst.score) {
                        arr2.remove(max-1);
                        add2 = true;
                    }
                }
            }
            else {
                arr2 = new ArrayList<>(max);
                add2 = true;
                neighbors.put(name2, arr2);
            }
            
            if (add1 || add2) {
                NeighborPair p = new NeighborPair(name1, name2, score);
                
                if (add1) {
                    arr1.add(p);
                    if (arr1.size() > 1) {
                        Collections.sort(arr1);
                    }
                }

                if (add2) {
                    arr2.add(p);
                    if (arr2.size() > 1) {
                        Collections.sort(arr2);
                    }
                }
            }            
        }
    }
    
    private class BestNeighbors {
        public int max = 6;
        
        public HashMap<String, ArrayList<Neighbor>> neighbors = new HashMap<>();
        
        public BestNeighbors() {
        }
        
        public BestNeighbors(int max) {
            this.max = max;
        }
        
        public void add(String name1, String name2, int score) {
            addHelper(name1, name2, score);
            addHelper(name2, name1, score);
        }
        
        private void addHelper(String target, String query, int score) {
            if (neighbors.containsKey(target)) {
                ArrayList<Neighbor> arr = neighbors.get(target);
                if (arr.size() < max) {
                    arr.add(new Neighbor(query, score));
                    Collections.sort(arr);
                }
                else {
                    Neighbor worst = arr.get(max-1);
                    if (score > worst.score) {
                        arr.set(max-1, new Neighbor(query, score));
                        Collections.sort(arr);
                    }
                }
            }
            else {
                ArrayList<Neighbor> arr = new ArrayList<>();
                arr.add(new Neighbor(query, score));
                neighbors.put(target, arr);
            }
        }
        
        public HashSet<String> getConnectedNeighbors(final String target, 
                HashSet<String> visited, HashSet<String> ignored) {
            
            visited.add(target);
            
            HashSet<String> connectedNeighbors = new HashSet<>();
            
            HashSet<String> pending = new HashSet<>();
            
            for (Neighbor n : neighbors.get(target)) {
                String name = n.name;
                if (!ignored.contains(name)) {
                    if (visited.contains(name)) {
                        connectedNeighbors.add(name);
                    }
                    else {
                        pending.add(name);
                    }
                }
            }
            
            while (!pending.isEmpty()) {
                Iterator<String> itr = pending.iterator();
                String name = itr.next();
                itr.remove();
                
                connectedNeighbors.add(name);
                visited.add(name);

                for (Neighbor n : neighbors.get(name)) {
                    String name2 = n.name;
                    if (!ignored.contains(name2)) {
                        if (visited.contains(name2)) {
                            connectedNeighbors.add(name2);
                        }
                        else {
                            pending.add(name2);
                        }
                    }
                }
            }
            
            return connectedNeighbors;
        }
    }
    
    private CONTAIN_STATUS getContained(PafRecord r) {
        return rnabloom.util.PafUtils.getContained(r, maxEdgeClip);
    }
    
    private CONTAIN_STATUS getContained(OverlapCoords r,
            int qMinStart, int qMaxEnd, int tMinStart, int tMaxEnd) {
        
        boolean qContained = r.qStart <= qMinStart + maxEdgeClip && qMaxEnd - r.qEnd <= maxEdgeClip;
        boolean tContained = r.tStart <= tMinStart + maxEdgeClip && tMaxEnd - r.tEnd <= maxEdgeClip;
        
        if (qContained && tContained) {
            int qLeftOver = (r.qStart - qMinStart) + (qMaxEnd - r.qEnd);
            int tLeftOver = (r.tStart - tMinStart) + (tMaxEnd - r.tEnd);
            if (qLeftOver < tLeftOver) {
                return CONTAIN_STATUS.QUERY;
            }
            else {
                return CONTAIN_STATUS.TARGET;
            }
        }
        else if (qContained) {
            return CONTAIN_STATUS.QUERY;
        }
        else if (tContained) {
            return CONTAIN_STATUS.TARGET;
        }
        
        return CONTAIN_STATUS.NEITHER;
    }

    private void findContained(ArrayDeque<Overlap> records, Set<String> contained, HashMap<String, Histogram> histogramMap) {        
        records.parallelStream().forEach(r -> {
            Histogram qHist = histogramMap.get(r.qName);
            Histogram tHist = histogramMap.get(r.tName);
            CONTAIN_STATUS c = getContained(r, qHist.minStart, qHist.maxEnd, tHist.minStart, tHist.maxEnd);
            switch (c) {
                case QUERY:
                    contained.add(r.qName);
                    break;
                case TARGET:
                    contained.add(r.tName);
                    break;
            }
        });
    }
        
    private void findContainedTargetOverlaps(String qName, Histogram qHist,
            ArrayDeque<TargetOverlap> records, Set<String> contained, HashMap<String, Histogram> histogramMap) {
        records.parallelStream().forEach(r -> {
            Histogram tHist = histogramMap.get(r.tName);
            CONTAIN_STATUS c = getContained(r, qHist.minStart, qHist.maxEnd, tHist.minStart, tHist.maxEnd);
            switch (c) {
                case TARGET:
                    // designate as "contained" only if the other is unique, i.e. not "contained"
                    if (!contained.contains(qName)) {
                        contained.add(r.tName);
                        // histogram of "contained" becomes irrelevant as it will not be evaluated
                        tHist.bars = null;
                    }
                    break;
                case QUERY:
                    if (!contained.contains(r.tName)) {
                        contained.add(qName);
                        qHist.bars = null;
                    }
                    break;
            }
        });
    }
    
    private void findContainedQueryOverlaps(String tName, Histogram tHist,
            Set<String> contained, HashMap<String, Histogram> histogramMap) {                
        tHist.pendingQueries.parallelStream().forEach(r -> {
            Histogram qHist = histogramMap.get(r.qName);
            CONTAIN_STATUS c = getContained(r, qHist.minStart, qHist.maxEnd, tHist.minStart, tHist.maxEnd);
            switch (c) {
                case QUERY:
                    // designate as "contained" only if the other is unique, i.e. not "contained"
                    if (!contained.contains(tName)) {
                        contained.add(r.qName);
                        // histogram of "contained" becomes irrelevant as it will not be evaluated
                        qHist.bars = null;
                    }
                    break;
                case TARGET:
                    if (!contained.contains(r.qName)) {
                        contained.add(tName);
                        tHist.bars = null;
                    }
                    break;
            }
        });
    }
    
    private void findContainedTargetOverlaps(String qName, Histogram qHist,
            ArrayDeque<TargetOverlap> records, Set<String> contained, 
            HashMap<String, Histogram> histogramMap,
            HashMap<String, PolyAInfo> polyAInfoMap) {
        
        PolyAInfo qInfo = polyAInfoMap.get(qName);
        
        records.parallelStream().forEach(r -> {
            Histogram tHist = histogramMap.get(r.tName);
            CONTAIN_STATUS c = getContained(r, qHist.minStart, qHist.maxEnd, tHist.minStart, tHist.maxEnd);
            switch (c) {
                case TARGET:
                    // designate as "contained" only if the other is unique, i.e. not "contained"
                    if (!contained.contains(qName)) {
                        PolyAInfo tInfo = polyAInfoMap.get(r.tName);
                        if (tInfo == null || isTargetPolyATContained(r, tInfo)) {
                            contained.add(r.tName);
                            // histogram of "contained" becomes irrelevant as it will not be evaluated
                            tHist.bars = null;
                        }
                    }
                    break;
                case QUERY:
                    if (!contained.contains(r.tName)) {
                        if (qInfo == null || isQueryPolyATContained(r, qInfo)) {
                            contained.add(qName);
                            qHist.bars = null;
                        }
                    }
                    break;
            }
        });
    }
    
    private void findContainedQueryOverlaps(String tName, Histogram tHist,
            Set<String> contained, HashMap<String, Histogram> histogramMap,
            HashMap<String, PolyAInfo> polyAInfoMap) {
        
        PolyAInfo tInfo = polyAInfoMap.get(tName);
        
        tHist.pendingQueries.parallelStream().forEach(r -> {
            Histogram qHist = histogramMap.get(r.qName);
            CONTAIN_STATUS c = getContained(r, qHist.minStart, qHist.maxEnd, tHist.minStart, tHist.maxEnd);
            switch (c) {
                case QUERY:
                    // designate as "contained" only if the other is unique, i.e. not "contained"
                    if (!contained.contains(tName)) {
                        PolyAInfo qInfo = polyAInfoMap.get(r.qName);
                        if (qInfo == null || isQueryPolyATContained(r, qInfo)) {
                            contained.add(r.qName);
                            // histogram of "contained" becomes irrelevant as it will not be evaluated
                            qHist.bars = null;
                        }
                    }
                    break;
                case TARGET:
                    if (!contained.contains(r.qName)) {
                        if (tInfo == null || isTargetPolyATContained(r, tInfo)) {
                            contained.add(tName);
                            tHist.bars = null;
                        }
                    }
                    break;
            }
        });
    }
    
    public long extractUniqueFromOverlaps(String outFastaPath) throws IOException {
        Timer timer = new Timer();
        
        Set<String> containedSet = Collections.synchronizedSet(new HashSet<>());
        HashMap<String, Histogram> histogramMap = new HashMap<>(100000);
        final int minSegmentLength = minOverlapMatches;
//        final int maxOverlapSizeDiff = maxIndelSize;
        
        PafReader reader = new PafReader(overlapPafInputStream);
        final boolean checkNumAltReads = minNumAltReads > 0;
        
        ArrayDeque<TargetOverlap> currRecords = new ArrayDeque<>();
        String currName = null;
        Histogram currHist = null;
        int currHistBinSize = -1;
        boolean currContained = false;
        
        for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if ((!stranded || !r.reverseComplemented) &&
                    !r.qName.equals(r.tName) &&
                    hasAlignment(r) &&
                    (hasLargeOverlap(r) || isContainmentPafRecord(r))) {
//                    (hasLargeOverlap(r) || isContainmentPafRecord(r)) &&
//                    hasGoodOverlap(r)) {
                
                //boolean hasAln = hasAlignment(r);
                //if ((!hasAln && hasSimilarSizedOverlap(r, maxOverlapSizeDiff)) || (hasAln && hasGoodAlignment(r))) {
                    if (!r.qName.equals(currName)) {
                        if (currName != null && currHist != null) {
                            if (!currRecords.isEmpty()) {
                                if (polyAInfoMap != null) {
                                    findContainedTargetOverlaps(currName, currHist, currRecords,
                                            containedSet, histogramMap, polyAInfoMap);
                                }
                                else {
                                    findContainedTargetOverlaps(currName, currHist, currRecords,
                                            containedSet, histogramMap);
                                }
                                currRecords.clear();
                            }

                            if (currHist.pendingQueries != null) {
                                if (polyAInfoMap != null) {
                                    findContainedQueryOverlaps(currName, currHist, containedSet,
                                        histogramMap, polyAInfoMap);
                                }
                                else {
                                    findContainedQueryOverlaps(currName, currHist, containedSet,
                                        histogramMap);
                                }
                                currHist.pendingQueries = null;
                            }
                        }

                        currName = r.qName; 

                        currHist = histogramMap.get(currName);
                        currHistBinSize = getHistogramBinSize(r.qLen);
                        if (currHist == null) {
                            currHist = new Histogram(r.qLen, r.qStart, r.qEnd, currHistBinSize);
                            histogramMap.put(currName, currHist);
                        }
                        currHist.seenAsQuery = true;

                        currContained = containedSet.contains(currName);
                    }

                    Histogram tHist = histogramMap.get(r.tName);
                    int tHistBinSize = getHistogramBinSize(r.tLen);
                    if (tHist == null) {
                        tHist = new Histogram(r.tLen, r.tStart, r.tEnd, tHistBinSize);
                        histogramMap.put(r.tName, tHist);
                    }

                    ArrayDeque<Interval>[] qtBlocks = getAlignedBlocks(r, maxIndelSize);
                    
                    // update coverage histogram for query
                    for (Interval qBlock : qtBlocks[0]) {
                        updateHistogram(currHist, qBlock.start, qBlock.end, currHistBinSize);
                    }
                    
                    // update coverage histogram for target
                    for (Interval tBlock : qtBlocks[1]) {
                        updateHistogram(tHist, tBlock.start, tBlock.end, tHistBinSize);
                    }
                    
                    //updateHistogram(currHist, r.qStart, r.qEnd, currHistBinSize);
                    //updateHistogram(tHist, r.tStart, r.tEnd, tHistBinSize);

                    if (qtBlocks[0].size() == 1 && qtBlocks[1].size() == 1 && // query and target overlap as a single-block alignment
                            !currContained && !containedSet.contains(r.tName)) { // query and target are not already designated as "contained"
                        // look for containment between query and target
                        if (tHist.seenAsQuery || isFullyCovered(tHist)) {
                            currRecords.add(pafToTargetOverlap(r));
                        }
                        else {
                            ArrayDeque<QueryOverlap> pending = tHist.pendingQueries;
                            if (pending == null) {
                                pending = new ArrayDeque<>();
                                tHist.pendingQueries = pending;
                            }
                            pending.add(pafToQueryOverlap(r));
                        }
                    }
                //}
            }
        }
        reader.close();
        
        if (!currRecords.isEmpty() && currHist != null) {
            if (polyAInfoMap != null) {
                findContainedTargetOverlaps(currName, currHist, currRecords,
                        containedSet, histogramMap, polyAInfoMap);
            }
            else {
                findContainedTargetOverlaps(currName, currHist, currRecords,
                        containedSet, histogramMap);
            }
            currRecords.clear();
        }
        
        for (HashMap.Entry<String, Histogram> e : histogramMap.entrySet()) {
            Histogram h = e.getValue();
            if (h.pendingQueries != null) {
                if (polyAInfoMap != null) {
                    findContainedQueryOverlaps(e.getKey(), h, containedSet,
                            histogramMap, polyAInfoMap);
                }
                else {
                    findContainedQueryOverlaps(e.getKey(), h, containedSet,
                            histogramMap);
                }
                h.pendingQueries = null;
            }
        }
        
        printMessage("Parsed " + NumberFormat.getInstance().format(reader.getNumRecords()) + " overlap records in " + timer.elapsedDHMS());
        
        // look for unique reads
        timer.start();
        Set<String> readsWithOverlap = histogramMap.keySet();
                
        FastaReader fr = new FastaReader(seqFastaPath);
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        long seqID = 0;
        long originalNumSeq = 0;
        long numMultiSeqs = 0;
        for (FastaRecord record = new FastaRecord(); fr.hasNext();) {
            ++originalNumSeq;
            fr.nextWithName(record);
            String seqName = record.name;
            if (!containedSet.contains(seqName)) {
                if (checkNumAltReads) {
                    if (readsWithOverlap.contains(seqName)) {
                        Histogram hist = histogramMap.get(seqName);
                        if (hist != null) {                            
                            ArrayDeque<Interval> spans = extractEffectiveIntervals(hist, getHistogramBinSize(hist.length), minNumAltReads, minSegmentLength);
                            if (!spans.isEmpty()) {
                                ++seqID;
                                int seqLen = record.seq.length();
                                if (spans.size() > 1) {
                                    ++numMultiSeqs;
                                    String namePrefix = seqID + "_p";
                                    int segmentID = 1;
                                    for (Interval i : spans) {
                                        String seg = record.seq.substring(i.start, i.end);
                                        if (!isLowComplexityLongWindowed(seg)) {
                                            String header = namePrefix + segmentID + " " + seqLen + ":" + i.start + "-" + i.end + " " + seqName;
                                            fw.write(header, seg);
                                            ++segmentID;
                                        }
                                    }
                                }
                                else {
                                    Interval i = spans.peekFirst();
                                    String seg = record.seq.substring(i.start, i.end);
                                    if (!isLowComplexityLongWindowed(seg)) {
                                        String header = Long.toString(seqID) + " " + seqLen + ":" + i.start + "-" + i.end  + " " + seqName;
                                        fw.write(header, seg);
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    fw.write(Long.toString(++seqID), record.seq);
                }
            }
        }
        fr.close();
        fw.close();
        
        printMessage("total reads:    " + NumberFormat.getInstance().format(originalNumSeq));
        printMessage(" - unique:      " + NumberFormat.getInstance().format(seqID) + "\t(" + convertToRoundedPercent(seqID/(float)originalNumSeq) + " %)");
        printMessage("   - multi-seg: " + NumberFormat.getInstance().format(numMultiSeqs));
        printMessage("Unique reads extracted in " + timer.elapsedDHMS());
        
        return seqID;
    }
    
    public class SeededCluster {
        final HashSet<String> seedNames = new HashSet<>();
        final HashSet<String> readNames = new HashSet<>();
        
        public void addSeedName(String name) {
            seedNames.add(name);
        }

        public void addSeedNames(Collection<String> names) {
            seedNames.addAll(names);
        }
        
        public void addReadName(String name) {
            readNames.add(name);
        }
        
        public int size() {
            return readNames.size();
        }
        
        public void merge(SeededCluster other) {
            this.seedNames.addAll(other.seedNames);
            this.readNames.addAll(other.readNames);
        }
    }
    
    private void addReadName(String readName, ArrayList<Neighbor> targets,
            HashMap<String, SeededCluster> seedNameClusterMap,
            HashMap<String, Integer> seedPairCounts) {
        
        String seedName = null;
        
        if (targets.size() == 1) {
            seedName = targets.iterator().next().name;
        }
        else {
            Collections.sort(targets);
            seedName = targets.get(0).name;
            String tuple = getMergedKey(seedName, targets.get(1).name);
            Integer count = seedPairCounts.get(tuple);
            if (count == null) {
                seedPairCounts.put(tuple, 1);
            }
            else {
                seedPairCounts.put(tuple, count+1);
            }
        }
        
        SeededCluster cluster = seedNameClusterMap.get(seedName);        
        if (cluster == null) {
            cluster = new SeededCluster();
            cluster.addSeedName(seedName);
            seedNameClusterMap.put(seedName, cluster);
        }
        
        cluster.addReadName(readName);
    }
    
    private void mergeClusters(HashMap<String, SeededCluster> seedNameClusterMap) {
        HashSet<String> visited = new HashSet<>();
        
        for (String seed : new ArrayDeque<>(seedNameClusterMap.keySet())) {
            if (!visited.contains(seed)) {
                SeededCluster mainCluster = seedNameClusterMap.get(seed);
                visited.add(seed);
                
                HashSet<String> pendingSeedNames = new HashSet<>(mainCluster.seedNames);
                pendingSeedNames.removeAll(visited);
                
                while (!pendingSeedNames.isEmpty()) {
                    Iterator<String> itr = pendingSeedNames.iterator();
                    String s = itr.next();
                    itr.remove();
                    
                    SeededCluster c = seedNameClusterMap.get(s);
                    visited.add(s);
                    
                    if (c == null) {
                        seedNameClusterMap.put(s, mainCluster);
                    } 
                    else if (c != mainCluster) {
                        HashSet<String> update = new HashSet<>(c.seedNames);
                        update.removeAll(visited);
                        pendingSeedNames.addAll(update);
                        mainCluster.merge(c);
                        
                        seedNameClusterMap.put(s, mainCluster);
                    }
                }
            }
        }
    }
    
    public void trimSplitByReadDepth(String targetFastaPath, String outFastaPath, int minSegmentLength) throws IOException {
        Timer timer = new Timer();
        
        HashMap<String, Histogram> histogramMap = new HashMap<>();
        
        PafReader reader = new PafReader(overlapPafInputStream);
        for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if ((!stranded || !r.reverseComplemented) &&
                    hasLargeOverlap(r) && //r.numMatch >= minOverlapMatches && 
                    (!hasAlignment(r) || hasGoodAlignment(r) || hasGoodMatches(r, minAlnId))) {                
                Histogram h = histogramMap.get(r.tName);
                int binSize = getHistogramBinSize(r.tLen);
                if (h == null) {
                    h = new Histogram(r.tLen, r.tStart, r.tEnd, binSize);
                    histogramMap.put(r.tName, h);
                }
                else {
                    updateHistogram(h, r.tStart, r.tEnd, binSize);
                }
            }
        }
        reader.close();
        
        printMessage("Parsed " + NumberFormat.getInstance().format(reader.getNumRecords()) + " mapping records in " + timer.elapsedDHMS());
        
        timer.start();
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        FastaReader fr = new FastaReader(targetFastaPath);
        FastaRecord record = new FastaRecord();
        int numSplitted = 0;
        int minDepth = minNumAltReads + 1;
        while (fr.hasNext()) {
            fr.nextWithName(record);
            
            Histogram hist = histogramMap.get(record.name);
            if (hist != null) {
                ArrayDeque<Interval> spans = extractEffectiveIntervals(hist, getHistogramBinSize(hist.length), minDepth, minSegmentLength);
                if (!spans.isEmpty()) {
                    if (spans.size() > 1) {
                        ++numSplitted;
                        String namePrefix = record.name + "_p";
                        int segmentID = 0;
                        for (Interval i : spans) {
                            if (i.end - i.start >= minSegmentLength) {
                                String seg = record.seq.substring(i.start, i.end);
                                String header = namePrefix + ++segmentID;
                                fw.write(header, seg);
                            }
                        }
                    }
                    else {
                        Interval i = spans.peekFirst();
                        if (i.end - i.start >= minSegmentLength) {
                            String seg = record.seq.substring(i.start, i.end);
                            String header = record.name;
                            fw.write(header, seg);
                        }
                    }
                }
            } 
        }
        fr.close();
        fw.close();
        
        printMessage("Splitted " + numSplitted + " reads.");
        printMessage("Wrote trimmed reads in " + timer.elapsedDHMS());
    }
    
    private String getMergedKey(String a, String b) {
        if (a.compareTo(b) > 0) {
            return b + " " + a;
        }
        else {
            return a + " " + b;
        }
    }
    
    private String[] splitMergedKey(String m) {
        return m.split(" ");
    }
    
    public void extractUniqueFromMapping(String outFastaPath) throws IOException {
        Timer timer = new Timer();
        
        HashMap<String, Histogram> histogramMap = new HashMap<>(100000);
        final int minSegmentLength = minOverlapMatches;
        final int maxOverlapSizeDiff = maxIndelSize;
        
        PafReader reader = new PafReader(overlapPafInputStream);
        final boolean checkNumAltReads = minNumAltReads > 0;
                
        for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if ((!stranded || !r.reverseComplemented) &&
                    (hasLargeOverlap(r) || isContainmentPafRecord(r)) &&
                    hasGoodOverlap(r)) {
                
                boolean hasAln = hasAlignment(r);
                if ((!hasAln && hasSimilarSizedOverlap(r, maxOverlapSizeDiff)) || (hasAln && hasGoodAlignment(r))) {
                    Histogram tHist = histogramMap.get(r.tName);
                    int tHistBinSize = getHistogramBinSize(r.tLen);
                    if (tHist == null) {
                        tHist = new Histogram(r.tLen, r.tStart, r.tEnd, tHistBinSize);
                        histogramMap.put(r.tName, tHist);
                    }
                    updateHistogram(tHist, r.tStart, r.tEnd, tHistBinSize);
                }
            }
        }
        reader.close();
        printMessage("Parsed " + NumberFormat.getInstance().format(reader.getNumRecords()) + " mapping records in " + timer.elapsedDHMS());
        
        // look for unique reads
        timer.start();
        Set<String> seqsWithOverlap = histogramMap.keySet();
                
        FastaReader fr = new FastaReader(seqFastaPath);
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        long seqID = 0;
        long originalNumSeq = 0;
        long numMultiSeqs = 0;
        int minDepth = minNumAltReads + 1;
        for (FastaRecord record = new FastaRecord(); fr.hasNext();) {
            ++originalNumSeq;
            fr.nextWithName(record);
            String seqName = record.name;
                if (checkNumAltReads) {
                    if (seqsWithOverlap.contains(seqName)) {
                        Histogram hist = histogramMap.get(seqName);
                        if (hist != null) {                            
                            ArrayDeque<Interval> spans = extractEffectiveIntervals(hist, getHistogramBinSize(hist.length), minDepth, minSegmentLength);
                            if (!spans.isEmpty()) {
                                ++seqID;
                                int seqLen = record.seq.length();
                                if (spans.size() > 1) {
                                    ++numMultiSeqs;
                                    String namePrefix = seqID + "_p";
                                    int segmentID = 1;
                                    for (Interval i : spans) {
                                        String seg = record.seq.substring(i.start, i.end);
                                        if (!isLowComplexityLongWindowed(seg)) {
                                            String header = namePrefix + segmentID + " " + seqLen + ":" + i.start + "-" + i.end + " " + seqName;
                                            fw.write(header, seg);
                                            ++segmentID;
                                        }
                                    }
                                }
                                else {
                                    Interval i = spans.peekFirst();
                                    String seg = record.seq.substring(i.start, i.end);
                                    if (!isLowComplexityLongWindowed(seg)) {
                                        String header = Long.toString(seqID) + " " + seqLen + ":" + i.start + "-" + i.end  + " " + seqName;
                                        fw.write(header, seg);
                                    }
                                }
                            }
                        }
                    }
                }
                else {
                    fw.write(Long.toString(++seqID), record.seq);
                }
        }
        fr.close();
        fw.close();
        
        printMessage("total reads:    " + NumberFormat.getInstance().format(originalNumSeq));
        printMessage(" - unique:      " + NumberFormat.getInstance().format(seqID) + "\t(" + convertToRoundedPercent(seqID/(float)originalNumSeq) + " %)");
        printMessage("   - multi-seg: " + NumberFormat.getInstance().format(numMultiSeqs));
        printMessage("Unique reads extracted in " + timer.elapsedDHMS());
    }
    
    public int[] extractClustersFromMapping(String outdir) throws IOException {
        // cluster sequences by mapping them to target seeds
        // merge clusters if sequences map to multiple seeds
        
        Timer timer = new Timer();
        
        HashMap<String, SeededCluster> seedNameClusterMap = new HashMap<>();
        HashMap<String, Integer> seedPairCounts = new HashMap<>();
        
        PafReader reader = new PafReader(overlapPafInputStream);
        String prevName = null;
        ArrayList<Neighbor> targets = new ArrayList<>();
        for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if (!r.qName.equals(prevName)) {
                if (!targets.isEmpty()) {
                    addReadName(prevName, targets, seedNameClusterMap, seedPairCounts);
                    targets.clear();
                }
                
                prevName = r.qName;
            }
            
            if ((!stranded || !r.reverseComplemented) &&
                    hasLargeOverlap(r) && r.numMatch >= minOverlapMatches &&
                    (!hasAlignment(r) || hasGoodAlignment(r) || hasGoodMatches(r, minAlnId))) {                
                targets.add(new Neighbor(r.tName, r.numMatch));
            }
        }
        reader.close();
        
        // process targets for the final query
        if (!targets.isEmpty()) {
            addReadName(prevName, targets, seedNameClusterMap, seedPairCounts);
            targets.clear();
        }
        
        printMessage("Parsed " + NumberFormat.getInstance().format(reader.getNumRecords()) + " mapping records in " + timer.elapsedDHMS());
        
        timer.start();
        
        // process counts
        int numSeedPairs = seedPairCounts.size();
        int numGoodSeedPairs = 0;
        int minPairingReads = this.minNumAltReads + 1;
        float minProportion = 0.01f;
        for (Entry<String, Integer> e : seedPairCounts.entrySet()) {
            String[] tuple = splitMergedKey(e.getKey());
            int count = e.getValue();
            if (count >= minPairingReads) {
                String a = tuple[0];
                String b = tuple[1];
                SeededCluster clusterA = seedNameClusterMap.get(a);
                SeededCluster clusterB = seedNameClusterMap.get(b);
                if (clusterA != null && clusterB != null) {
                    int sizeA = clusterA.size();
                    int sizeB = clusterB.size();
                                
                    if (minPairingReads >= minProportion * Math.max(sizeA, sizeB)) {
                        ++numGoodSeedPairs;
                        clusterA.addSeedName(b);
                        clusterB.addSeedName(a);
                    }
                }
            }
        }
        
        printMessage("Kept " + numGoodSeedPairs + " of " + numSeedPairs +
                " (" + convertToRoundedPercent(numGoodSeedPairs/(float)numSeedPairs) + " %) seed pairs.");
        
        // merge overlapping clusters
        mergeClusters(seedNameClusterMap);
        
        // assign cluster IDs
        int numClusters = 0;
        int largestClusterSize = 0;
        int largestClusterID = -1;
        HashMap<String, Integer> readClusterIDs = new HashMap<>();
        String clusterDescPath = outdir + File.separator + "cluster_desc.txt";
        FileWriter writer = new FileWriter(clusterDescPath);
        writer.write("cid\tsize\tnum_seeds\tseeds\n");
        for (SeededCluster c : seedNameClusterMap.values()) {
            int size = c.size();
            if (size > 0) {       
                ++numClusters;

                for (String n : c.readNames) {
                    readClusterIDs.put(n, numClusters);
                }
                
                if (size > largestClusterSize) {
                    largestClusterSize = size;
                    largestClusterID = numClusters;
                }
                
                writer.write(numClusters + "\t" + c.readNames.size() + "\t" + c.seedNames.size() + "\t" + c.seedNames.toString() + "\n");
                
                // clear set because seeds in a merged cluster point to the same set
                c.readNames.clear();
            }
        }
        writer.close();
                
        printMessage("Cluster IDs assigned in " + timer.elapsedDHMS());
        printMessage("\t- clusters:\t" + numClusters);
        printMessage("\t- largest:\t#" + largestClusterID + " (" + largestClusterSize + " reads)");
        
        // create cluster directories
        // cluster "0" is for orphans
        for (int i=0; i<=numClusters; ++i) {
            String path = outdir + File.separator + i;
            new File(path).mkdirs();
        }
        
        System.gc();
        
        // write fasta files
        timer.start();
        FastaReader fr = new FastaReader(seqFastaPath);
        ArrayList<ArrayDeque<CompressedFastaRecord>> clusterRecords = new ArrayList<>(numClusters+1);
        clusterRecords.set(0, new ArrayDeque<>());
        while (fr.hasNext()) {
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            String seq = nameSeq[1];
            
            Integer cid = readClusterIDs.get(name);
            if (cid != null) {
                ArrayDeque<CompressedFastaRecord> fastaBuffer = clusterRecords.get(cid);
                if (fastaBuffer == null) {
                    fastaBuffer = new ArrayDeque<>();
                    clusterRecords.set(cid, fastaBuffer);
                }
                fastaBuffer.add(new CompressedFastaRecord(name, seq));
            } 
            else {
                clusterRecords.get(0).add(new CompressedFastaRecord(name, seq));
            }
        }
        fr.close();
        
        int[] counts = writeClusterFastaFiles(clusterRecords, outdir);
        System.gc();
        
        printMessage("\t- orphans:\t" + counts[0]);
        printMessage("Clustered reads written in " + timer.elapsedDHMS());
        printMessage("Clustering completed in " + timer.totalElapsedDHMS());
        
        return counts;
    }
    
    private int[] writeClusterFastaFiles(ArrayList<ArrayDeque<CompressedFastaRecord>> clusterRecords, String outdir) throws IOException {
        int size = clusterRecords.size();
        int[] counts = new int[size];
        for (int cid=0; cid<size; ++cid) {
            String filePath = outdir + File.separator + cid + File.separator + cid + FASTA_EXT + GZIP_EXT;
            FastaWriter fw = new FastaWriter(filePath, false);
            ArrayDeque<CompressedFastaRecord> records = clusterRecords.get(cid);
            for (CompressedFastaRecord f : records) {
                fw.write(f.name, f.seqbits.toString());
            }
            fw.close();
            counts[cid] = records.size();
        }
        return counts;
    }
    
    public int[] extractClustersFromOverlaps(String outdir, int maxMergedClusterSize) throws IOException {
        Timer timer = new Timer();
        PafReader reader = new PafReader(overlapPafInputStream);
        
        final boolean checkNumAltReads = minNumAltReads > 0;
        
        HashMap<String, short[]> readHistogramMap = new HashMap<>();
        final int histBinSize = 25;
        final int maxBestNeighbors = 2;
        BestNeighborPairs bestNeighbors = new BestNeighborPairs(maxBestNeighbors);
//        BestNeighbor bestNeighbors = new BestNeighbor();
        
        for (PafRecord r = new PafRecord(); reader.hasNext();) {            
            reader.next(r);
            
            if ((!stranded || !r.reverseComplemented) &&
                    hasLargeOverlap(r) && hasGoodOverlap(r) &&
                    !r.qName.equals(r.tName)) {
                if (checkNumAltReads) {
                    short[] hist = readHistogramMap.get(r.qName);
                    if (hist == null) {
                        hist = initShortHistogram(r.qLen, histBinSize);
                        readHistogramMap.put(r.qName, hist);
                    }
                    updateHistogram(hist, histBinSize, r.qLen, r.qStart, r.qEnd);
                    
                    hist = readHistogramMap.get(r.tName);
                    if (hist == null) {
                        hist = initShortHistogram(r.tLen, histBinSize);
                        readHistogramMap.put(r.tName, hist);
                    }
                    updateHistogram(hist, histBinSize, r.tLen, r.tStart, r.tEnd);
                }
                
//                if (isDovetailPafRecord(r) || isContainmentPafRecord(r)) {
//                    bestNeighbors.add(r.qName, r.tName, r.numMatch);
//                }
                bestNeighbors.add(r.qName, r.tName, r.numMatch);
            }
        }
        reader.close();
        
        printMessage("Parsed " + NumberFormat.getInstance().format(reader.getNumRecords()) + " overlap records in " + timer.elapsedDHMS());
        
        // identify multi-segment reads and extract effective intervals
        final int minHistSegLen = minOverlapMatches/histBinSize;
        HashSet<String> multiSegmentSeqs = new HashSet<>();
        ConcurrentHashMap<String, ArrayDeque<Interval>> readSpansMap = new ConcurrentHashMap<>(bestNeighbors.neighbors.size());
        
        if (checkNumAltReads) {
            timer.start();
            Set<String> syncSet = Collections.synchronizedSet(multiSegmentSeqs);

            readHistogramMap.entrySet().parallelStream().forEach(
                e -> {
                    String name = e.getKey();
                    short[] hist = e.getValue();
                    if (hist != null) {
                        ArrayDeque<Interval> spans = extractEffectiveIntervals(hist, histBinSize, minNumAltReads, minHistSegLen);
                        readSpansMap.put(name, spans);
                        if (spans.size() > 1) {
                            syncSet.add(name);
                        }
                    }
                }
            );

            printMessage("Effective intervals identified in " + timer.elapsedDHMS());
            printMessage("\t- multi-segs:\t" + multiSegmentSeqs.size());
        }
        
        timer.start();
        ArrayList<NeighborPair> neighborPairsList = new ArrayList<>((bestNeighbors.neighbors.size() - multiSegmentSeqs.size()) * 2);
        ArrayList<NeighborPair> multiSegNeighborPairsList = new ArrayList<>(multiSegmentSeqs.size() * 2);
        /**@TODO split list into tiers to avoid list size reaching max Integer value*/
        for (ArrayList<NeighborPair> a : bestNeighbors.neighbors.values()) {
            for (NeighborPair p : a) {
                if (!multiSegmentSeqs.contains(p.name1) && !multiSegmentSeqs.contains(p.name2)) {
                    neighborPairsList.add(p);
                }
                else {
                    multiSegNeighborPairsList.add(p);
                }
            }
        }
        
        Collections.sort(neighborPairsList);
        Collections.sort(multiSegNeighborPairsList);
        
        printMessage("Overlaps sorted in " + timer.elapsedDHMS());
        
        // form clusters
        timer.start();
        ReadClusters3 clusters = new ReadClusters3(maxMergedClusterSize);
        
        NeighborPair last = null; // sorted list may contain duplicates
        for (NeighborPair p : neighborPairsList) {
            if (p != last) {
                clusters.add(p.name1, p.name2);
            }
            last = p;
        }
        
        last = null; // sorted list may contain duplicates
        for (NeighborPair p : multiSegNeighborPairsList) {
            if (p != last) {
                clusters.add(p.name1, p.name2);
            }
            last = p;
        }
        
        int numClusters = clusters.size();
        
        printMessage("Clusters formed in " + timer.elapsedDHMS());
        printMessage("\t- clusters:\t" + numClusters);
        
        // assign cluster IDs
        timer.start();
        clusters.assignIDs();
        int largestClusterID = clusters.getLargestClusterID();
        int largestClusterSize= clusters.getLargestClusterSize();
        printMessage("Cluster IDs assigned in " + timer.elapsedDHMS());
        printMessage("\t- largest:\t#" + largestClusterID + " (" + largestClusterSize + " reads)");

        // extract effective regions for each read
        timer.start();
        FastaReader fr = new FastaReader(seqFastaPath);
        int[] counts = new int[numClusters];
        int numOrphans = 0;
        final int maxBufferSize = 1000000;
        int numReadsInBuffer = 0;
        ArrayDeque<CompressedFastaRecord> orphanRecords = new ArrayDeque<>();
        HashMap<Integer, ArrayDeque<CompressedFastaRecord>> clusterRecords = new HashMap<>(Math.min(numClusters, maxBufferSize));
        while (fr.hasNext()) {
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            String seq = nameSeq[1];
            
            int cid = clusters.getID(name);
            ArrayDeque<CompressedFastaRecord> fastaBuffer = null;
            if (cid > 0) {
                counts[cid-1] += 1;
                fastaBuffer = clusterRecords.get(cid);
                if (fastaBuffer == null) {
                    fastaBuffer = new ArrayDeque<>();
                    clusterRecords.put(cid, fastaBuffer);
                }
                ++numReadsInBuffer;
            } 
            else {
                ++numOrphans;
                if (!checkNumAltReads) {
                    // number of alternate reads requires is 0
                    // sequences with read depth of 1 (itself) will be kept
                    orphanRecords.add(new CompressedFastaRecord(name, seq));
                }
                continue;
            }
            
            ArrayDeque<Interval> spans = readSpansMap.get(name);
            if (checkNumAltReads && spans != null) {
                int numSpans = spans.size();
                int seqLen = seq.length();
                String seqLenStr = Integer.toString(seqLen);
                if (numSpans == 1) {
                    Interval span = spans.peekFirst();
                    int start = Math.max(0, span.start);
                    int end = Math.min(span.end, seqLen);
                    String header = name + "_t " + seqLenStr + ":" + start + "-" + end;
                    if (start == 0 && end == seqLen) {
                        fastaBuffer.add(new CompressedFastaRecord(header, seq));
                    }
                    else {
                        fastaBuffer.add(new CompressedFastaRecord(header, seq.substring(start, end)));
                    }
                }
                else {
                    int i = 1;
                    for (Interval span : spans) {
                        int start = Math.max(0, span.start);
                        int end = Math.min(span.end, seqLen);
                        String header = name + "_p" + i++ + " " + seqLenStr + ":" + start + "-" + end;
                        fastaBuffer.add(new CompressedFastaRecord(header, seq.substring(start, end)));
                    }
                }
            }
            else {
                fastaBuffer.add(new CompressedFastaRecord(name, seq));
            }
            
            if (numReadsInBuffer >= maxBufferSize) {
                emptyClusterFastaBuffer(clusterRecords, outdir, true);
                numReadsInBuffer = 0;
            }
        }
        fr.close();
        
        emptyClusterFastaBuffer(clusterRecords, outdir, true);
        
        printMessage("\t- orphans:\t" + numOrphans);
        if (!checkNumAltReads) {
            emptyOrphanRecords(orphanRecords, outdir, false);
        }
        
        printMessage("Clustered reads written in " + timer.elapsedDHMS());
        printMessage("Clustering completed in " + timer.totalElapsedDHMS());
        System.gc();
        
        //System.out.println("before: " + NumberFormat.getInstance().format(originalNumSeq) + "\tafter: " + NumberFormat.getInstance().format(seqID));
        return counts;
    }
    
    private void emptyClusterFastaBuffer(HashMap<Integer, ArrayDeque<CompressedFastaRecord>> clusterRecords, 
            String outdir, boolean append) throws IOException {
        if (!clusterRecords.isEmpty()) {
            for (Map.Entry<Integer, ArrayDeque<CompressedFastaRecord>> e : clusterRecords.entrySet()) {
                int cid = e.getKey();
                String filePath = outdir + File.separator + cid + File.separator + cid + FASTA_EXT + GZIP_EXT;
                FastaWriter fw = new FastaWriter(filePath, append);
                for (CompressedFastaRecord f : e.getValue()) {
                    fw.write(f.name, f.seqbits.toString());
                }
                fw.close();
            }
            clusterRecords.clear();
        }
    }
    
    private void emptyOrphanRecords(ArrayDeque<CompressedFastaRecord> records, 
            String outdir, boolean append) throws IOException {
        String filePath = outdir + File.separator + "0" + File.separator + "0" + FASTA_EXT + GZIP_EXT;
        FastaWriter fw = new FastaWriter(filePath, append);
        for (CompressedFastaRecord f : records) {
            fw.write(f.name, f.seqbits.toString());
        }
        fw.close();
        records.clear();
    }
        
    private void removeVertex(String name) {
        removeVertexStranded(name + '+');
        removeVertexStranded(name + '-');
    }
    
    private void removeVertexStranded(String strandedName) {
        if (graph.containsVertex(strandedName)) {
            graph.removeAllEdges(graph.edgesOf(strandedName));
            graph.removeVertex(strandedName);
        }
    }
    
    private void addVertex(String name) {
        addVertexStranded(name + '+');
        addVertexStranded(name + '-');
    }
    
    private void addVertexStranded(String strandedName) {
        if (!graph.containsVertex(strandedName)) {
            graph.addVertex(strandedName);
        }
    }
    
    private void addEdges(PafRecord r) {
        if (r.reverseComplemented) {
            if (r.qEnd >= r.qLen - maxEdgeClip && r.tEnd >= r.tLen - maxEdgeClip) {
                addVertex(r.qName);
                addVertex(r.tName);
                graph.addEdge(r.tName+"+", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                graph.addEdge(r.qName+"+", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
            }
            else if (r.tStart <= maxEdgeClip && r.qStart <= maxEdgeClip) {
                addVertex(r.qName);
                addVertex(r.tName);
                graph.addEdge(r.qName+"-", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                graph.addEdge(r.tName+"-", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
            }
        }
        else {
            if (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) {
                addVertex(r.qName);
                addVertex(r.tName);
                graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                graph.addEdge(r.tName+"-", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
            }
            else if (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip) {
                addVertex(r.qName);
                addVertex(r.tName);
                graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                graph.addEdge(r.qName+"-", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
            }
        }
    }
    
    private void addForwardEdge(PafRecord r) {
        if (!r.reverseComplemented) {
            if (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) {
                addVertexStranded(r.qName+"+");
                addVertexStranded(r.tName+"+");
                graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
            }
            else if (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip) {
                addVertexStranded(r.qName+"+");
                addVertexStranded(r.tName+"+");
                graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
            }
        }
    }
    
    private void addEdges(StrandedOverlap r) {
//        if (r.reverseComplemented) {
//            if (r.qEnd >= r.qLen - maxEdgeClip && r.tEnd >= r.tLen - maxEdgeClip) {
//                addVertex(r.qName);
//                addVertex(r.tName);
//                graph.addEdge(r.tName+"+", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
//                graph.addEdge(r.qName+"+", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
//            }
//            else if (r.tStart <= maxEdgeClip && r.qStart <= maxEdgeClip) {
//                addVertex(r.qName);
//                addVertex(r.tName);
//                graph.addEdge(r.qName+"-", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
//                graph.addEdge(r.tName+"-", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
//            }
//        }
//        else {
//            if (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) {
//                addVertex(r.qName);
//                addVertex(r.tName);
//                graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
//                graph.addEdge(r.tName+"-", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
//            }
//            else if (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip) {
//                addVertex(r.qName);
//                addVertex(r.tName);
//                graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
//                graph.addEdge(r.qName+"-", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
//            }
//        }
        
        int qHead = r.qStart;
        int qTail = r.qLen - r.qEnd;
        int tHead = r.tStart;
        int tTail = r.tLen - r.tEnd;
        
        if (r.reverseComplemented) {
            if (qHead > tTail && tHead > qTail && tTail <= maxEdgeClip && qTail <= maxEdgeClip) {
                addVertex(r.qName);
                addVertex(r.tName);
                graph.addEdge(r.tName+"+", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                graph.addEdge(r.qName+"+", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
            }
            else if (qTail > tHead && tTail > qHead && qHead <= maxEdgeClip && tHead <= maxEdgeClip) {
                addVertex(r.qName);
                addVertex(r.tName);
                graph.addEdge(r.qName+"-", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                graph.addEdge(r.tName+"-", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
            }
        }
        else {
            if (qHead > tHead && qTail < tTail && tHead <= maxEdgeClip && qTail <= maxEdgeClip) {
                addVertex(r.qName);
                addVertex(r.tName);
                graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                graph.addEdge(r.tName+"-", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
            }
            else if (tHead > qHead && tTail < qTail && qHead <= maxEdgeClip && tTail <= maxEdgeClip) {
                addVertex(r.qName);
                addVertex(r.tName);
                graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                graph.addEdge(r.qName+"-", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
            }
        }
    }
    
    private void addEdges(StrandedOverlap r, PolyAInfo qInfo, PolyAInfo tInfo) {        
        int qHead = r.qStart;
        int qTail = r.qLen - r.qEnd;
        int tHead = r.tStart;
        int tTail = r.tLen - r.tEnd;
        
        if (r.reverseComplemented) {
            if (qHead > tTail && tHead > qTail && tTail <= maxEdgeClip && qTail <= maxEdgeClip) {
                // q+ -> t-
                boolean qPolyA = qInfo != null && qInfo.polyATail != null;
                boolean tPolyA = tInfo != null && tInfo.polyATail != null;
                if ((qInfo == null && tInfo == null) ||
                        (!qPolyA && !tPolyA) ||
                        (qPolyA && r.qEnd >= qInfo.polyATail.start + minPolyALength) ||
                        (tPolyA && r.tEnd >= tInfo.polyATail.start + minPolyALength)) {
                    addVertex(r.qName);
                    addVertex(r.tName);
                    graph.addEdge(r.tName+"+", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                    graph.addEdge(r.qName+"+", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                }
            }
            else if (qTail > tHead && tTail > qHead && qHead <= maxEdgeClip && tHead <= maxEdgeClip) {
                // q- -> t+
                boolean qPolyT = qInfo != null && qInfo.polyTHead != null;
                boolean tPolyT = tInfo != null && tInfo.polyTHead != null;
                if ((qInfo == null && tInfo == null) ||
                        (!qPolyT && !tPolyT) ||
                        (qPolyT && r.qStart <= qInfo.polyTHead.end - minPolyALength) || 
                        (tPolyT && r.tStart <= tInfo.polyTHead.end - minPolyALength)) {
                    addVertex(r.qName);
                    addVertex(r.tName);
                    graph.addEdge(r.qName+"-", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                    graph.addEdge(r.tName+"-", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                }
            }
        }
        else {
            if (qHead > tHead && qTail < tTail && tHead <= maxEdgeClip && qTail <= maxEdgeClip) {
                // q+ -> t+
                boolean qPolyA = qInfo != null && qInfo.polyATail != null;
                boolean tPolyT = tInfo != null && tInfo.polyTHead != null;
                if ((qInfo == null && tInfo == null) ||
                        (!qPolyA && !tPolyT) ||
                        (qPolyA && r.qEnd >= qInfo.polyATail.start + minPolyALength) ||
                        (tPolyT && r.tStart <= tInfo.polyTHead.end - minPolyALength)) {
                    addVertex(r.qName);
                    addVertex(r.tName);
                    graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                    graph.addEdge(r.tName+"-", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                }
            }
            else if (tHead > qHead && tTail < qTail && qHead <= maxEdgeClip && tTail <= maxEdgeClip) {
                // t+ -> q+
                boolean tPolyA = tInfo != null && tInfo.polyATail != null;
                boolean qPolyT = qInfo != null && qInfo.polyTHead != null;
                                
                if ((qInfo == null && tInfo == null) ||
                        (!tPolyA && !qPolyT) ||
                        (tPolyA && r.tEnd >= tInfo.polyATail.start + minPolyALength) ||
                        (qPolyT && r.qStart <= qInfo.polyTHead.end - minPolyALength)) {
                    addVertex(r.qName);
                    addVertex(r.tName);
                    graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                    graph.addEdge(r.qName+"-", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                }
            }
        }
    }
    
    private void addForwardEdge(Overlap r) {
//        if (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) {
//            addVertexStranded(r.qName+"+");
//            addVertexStranded(r.tName+"+");
//            graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
//        }
//        else if (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip) {
//            addVertexStranded(r.qName+"+");
//            addVertexStranded(r.tName+"+");
//            graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
//        }

        int qHead = r.qStart;
        int qTail = r.qLen - r.qEnd;
        int tHead = r.tStart;
        int tTail = r.tLen - r.tEnd;
        
        if (qHead > tHead && qTail < tTail && tHead <= maxEdgeClip && qTail <= maxEdgeClip) {
            addVertex(r.qName);
            addVertex(r.tName);
            graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
            graph.addEdge(r.tName+"-", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
        }
        else if (tHead > qHead && tTail < qTail && qHead <= maxEdgeClip && tTail <= maxEdgeClip) {
            addVertex(r.qName);
            addVertex(r.tName);
            graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
            graph.addEdge(r.qName+"-", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
        }
    }
    
    private void printCounts(HashMap<String, Float> counts) {
        for (Entry<String, Float> entry : counts.entrySet()) {
            System.out.println(entry.getKey() + '\t' + entry.getValue());
        }
    }
    
    private void printGraph() {
        for (OverlapEdge e : graph.edgeSet()) {
            System.out.println(graph.getEdgeSource(e) + " -> " + graph.getEdgeTarget(e) +
                    " [o=" + (e.sinkEnd-e.sinkStart + e.sourceEnd-e.sourceStart)/2 + "]");
        }
    }
    
    private void writeGraph(String path) throws IOException {
        Writer writer = getTextFileWriter(path, false);
        for (OverlapEdge e : graph.edgeSet()) {
            writer.write(graph.getEdgeSource(e) + " -> " + graph.getEdgeTarget(e) +
                    " [o=" + (e.sinkEnd-e.sinkStart + e.sourceEnd-e.sourceStart)/2 + "]\n");
        }
        writer.close();
    }
        
    private class TargetOverlapComparator implements Comparator<Overlap> {
        @Override
        public int compare(Overlap o1, Overlap o2) {
            return o1.tName.compareTo(o2.tName);
        }
    }
    
    private class OverlapSizeComparator implements Comparator<Overlap> {
        @Override
        public int compare(Overlap o1, Overlap o2) {
            // sort by overlap size in descending order
            return (o2.qEnd - o2.qStart) - (o1.qEnd - o1.qStart);
        }
    }
    
    private class OverlapQueryStartComparator implements Comparator<Overlap> {
        @Override
        public int compare(Overlap o1, Overlap o2) {
            // sort by qStart in ascending order
            return o1.qStart - o2.qStart;
        }
    }
    
    private void removeReverseComplementArtifacts(
            ArrayList<StrandedOverlap> overlaps,
            HashSet<String> containedSet) {
        
        int numOverlaps = overlaps.size();

        StrandedOverlap o1 = overlaps.get(0);
        for (int i=1; i<numOverlaps; ++i) {
            StrandedOverlap o2 = overlaps.get(i);
            String qName = o2.qName;
            String tName = o2.tName;
            if (o1.reverseComplemented != o2.reverseComplemented &&
                    o1.tName.equals(tName) && !containedSet.contains(tName)) {
                Interval q1 = new Interval(o1.qStart, o1.qEnd);
                Interval q2 = new Interval(o2.qStart, o2.qEnd);
                Interval t1 = new Interval(o1.tStart, o1.tEnd);
                Interval t2 = new Interval(o2.tStart, o2.tEnd);

                if (isContained(q1, q2)) {
                    containedSet.add(tName);
                }
                else if (isContained(t1, t2)) {
                    containedSet.add(qName);
                    break;
                }
                else {
                    int tLen = o1.tLen;
                    int qLen = o1.qLen;
                    int tUnmapped = tLen;
                    int qUnmapped = qLen;

                    Interval qMerged = merge(q1, q2);
                    if (qMerged == null) {
                        qUnmapped -= q1.end - q1.start;
                        qUnmapped -= q2.end - q2.start;
                    }
                    else {
                        qUnmapped -= qMerged.end - qMerged.start;
                    }

                    Interval tMerged = merge(t1, t2);
                    if (tMerged == null) {
                        tUnmapped -= t1.end - t1.start;
                        tUnmapped -= t2.end - t2.start;
                    }
                    else {
                        tUnmapped -= tMerged.end - tMerged.start;
                    }

                    if (tUnmapped > qUnmapped) {
                        int qMin = Math.min(q1.start, q2.start);
                        int qMax = Math.max(q1.end, q2.end);
                        if (qMin <= maxEdgeClip && qMax + maxEdgeClip >= qLen) {
                            containedSet.add(qName);
                            break;
                        }
                    }
                    else {
                        int tMin = Math.min(t1.start, t2.start);
                        int tMax = Math.max(t1.end, t2.end);
                        if (tMin <= maxEdgeClip && tMax + maxEdgeClip >= tLen) {
                            containedSet.add(tName);
                        }
                    }
                }
            }
            o1 = o2;
        }
    }
    
    private HashSet<String> populateGraphFromOverlaps() throws IOException {
        Timer timer = new Timer();
        
        HashSet<String> containedSet = new HashSet<>();
        PafReader reader = new PafReader(overlapPafInputStream);
        TargetOverlapComparator overlapComparator = new TargetOverlapComparator();
        
        // look for containment and overlaps
        if (stranded) {
            ArrayDeque<Overlap> overlaps = new ArrayDeque<>();
            ArrayDeque<Overlap> containedDoveTailOverlaps = new ArrayDeque<>();
            
            String prevName = null;
            boolean prevContained = false;
            ArrayDeque<Overlap> pendingOverlaps = new ArrayDeque<>();
            ArrayList<StrandedOverlap> nonStandardOverlaps = new ArrayList<>();
            
            for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
                reader.next(r);                
                if (r.reverseComplemented) {
                    if (this.cutRevCompArtifact && hasLargeOverlap(r) && hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
                        nonStandardOverlaps.add(pafToStrandedOverlap(r));
                    }
                }
                else if (!r.qName.equals(r.tName)) {
                    if ((hasLargeOverlap(r) || isContainmentPafRecord(r)) && hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
                        if (prevName == null) {
                            prevName = r.qName;
                        }
                        else if (!prevName.equals(r.qName)) {
                            if (!pendingOverlaps.isEmpty()) {
                                if (!prevContained) {
                                    overlaps.addAll(pendingOverlaps);
                                }
                                pendingOverlaps.clear();
                            }
                            
                            if (this.cutRevCompArtifact && !nonStandardOverlaps.isEmpty()) {
                                if (!prevContained) {
                                    // remove artifacts based on non standard overlaps
                                    nonStandardOverlaps.sort(overlapComparator);
                                    removeReverseComplementArtifacts(nonStandardOverlaps, containedSet);
                                }
                                nonStandardOverlaps.clear();
                            }
                            
                            prevName = r.qName;
                            prevContained = containedSet.contains(prevName);
                        }
                                                
                        if (!prevContained && !containedSet.contains(r.tName)) {
                            CONTAIN_STATUS c = getContained(r);
                            boolean isDoveTail = isForwardDovetailPafRecord(r);
                            switch (c) {
                                case QUERY:
                                    if (isDoveTail) {
                                        containedDoveTailOverlaps.add(pafToOverlap(r));
                                    }
                                    else {
                                        containedSet.add(r.qName);
                                        prevContained = true;
                                    }
                                    break;
                                case TARGET:
                                    if (isDoveTail) {
                                        containedDoveTailOverlaps.add(pafToOverlap(r));
                                    }
                                    else {
                                        containedSet.add(r.tName);
                                    }
                                    break;
                                case BOTH:
                                    if (isDoveTail) {
                                        containedDoveTailOverlaps.add(pafToOverlap(r));
                                    }
                                    else {
                                        if (r.qLen > r.tLen) {
                                            containedSet.add(r.tName);
                                        }
                                        else if (r.tLen > r.qLen) {
                                            containedSet.add(r.qName);
                                            prevContained = true;
                                        }
                                        else {
                                            if (r.qName.compareTo(r.tName) > 0) {
                                                containedSet.add(r.tName);
                                            }
                                            else {
                                                containedSet.add(r.qName);
                                                prevContained = true;
                                            }
                                        }
                                    }
                                    break;
                                case NEITHER:
                                    if (isDoveTail) {
                                        pendingOverlaps.add(pafToOverlap(r));
                                    }
                                    else {
                                        nonStandardOverlaps.add(pafToStrandedOverlap(r));
                                    }
                                    break;
                            }
                        }
                    }
                }
            }
            reader.close();

            if (prevName != null && !pendingOverlaps.isEmpty()) {
                if (!prevContained) {
                    overlaps.addAll(pendingOverlaps);
                }
                pendingOverlaps.clear();
            }
            
            if (this.cutRevCompArtifact && !nonStandardOverlaps.isEmpty()) {
                if (!prevContained) {
                    // remove artifacts based on non standard overlaps
                    nonStandardOverlaps.sort(overlapComparator);
                    removeReverseComplementArtifacts(nonStandardOverlaps, containedSet);
                }
                nonStandardOverlaps.clear();
            }
            
            for (Overlap r : overlaps) {
                if (!containedSet.contains(r.qName) &&
                    !containedSet.contains(r.tName)) {
                    addForwardEdge(r);
                }
            }
            
            for (Overlap r : containedDoveTailOverlaps) {
                String tName = r.tName;
                String qName = r.qName;
                if (!containedSet.contains(qName) &&
                    !containedSet.contains(tName)) {
                    qName = qName + '+';
                    tName = tName + '+';
                    if ((!graph.containsVertex(qName) || (graph.inDegreeOf(qName) == 0 || graph.outDegreeOf(qName) == 0)) &&
                        (!graph.containsVertex(tName) || (graph.inDegreeOf(tName) == 0 || graph.outDegreeOf(tName) == 0))) {
                        addForwardEdge(r);
                    }
                }
            }
        }
        else {
            ArrayDeque<StrandedOverlap> overlaps = new ArrayDeque<>();
            ArrayDeque<StrandedOverlap> containedDoveTailOverlaps = new ArrayDeque<>();
            ArrayList<StrandedOverlap> nonStandardOverlaps = new ArrayList<>();
            
            String prevName = null;
            boolean prevContained = false;
            ArrayDeque<StrandedOverlap> pendingOverlaps = new ArrayDeque<>();
            
            if (polyAInfoMap != null) {
                for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
                    reader.next(r);
                    if (!r.qName.equals(r.tName)) {
                        if ((hasLargeOverlap(r) || isContainmentPafRecord(r)) && hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {                        
                            if (prevName == null) {
                                prevName = r.qName;
                            }
                            else if (!prevName.equals(r.qName)) {
                                if (!pendingOverlaps.isEmpty()) { 
                                    if (!prevContained) {
                                        overlaps.addAll(pendingOverlaps);
                                    }
                                    pendingOverlaps.clear();
                                }

                                if (this.cutRevCompArtifact && !nonStandardOverlaps.isEmpty()) {
                                    if (!prevContained) {
                                        // remove artifacts based on non standard overlaps
                                        nonStandardOverlaps.sort(overlapComparator);
                                        removeReverseComplementArtifacts(nonStandardOverlaps, containedSet);
                                    }
                                    nonStandardOverlaps.clear();
                                }

                                prevName = r.qName;
                                prevContained = containedSet.contains(prevName);
                            }

                            if (!prevContained && !containedSet.contains(r.tName)) {
                                CONTAIN_STATUS c = getContained(r);
                                boolean isDoveTail = isDovetailPafRecord(r);
                                switch (c) {
                                    case QUERY:
                                        if (isDoveTail) {
                                            containedDoveTailOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        else {
                                            PolyAInfo qInfo = polyAInfoMap.get(r.qName);
                                            if (qInfo == null || isQueryPolyATContained(r, qInfo)) {
                                                containedSet.add(r.qName);
                                                prevContained = true;
                                            }
                                        }
                                        break;
                                    case TARGET:
                                        if (isDoveTail) {
                                            containedDoveTailOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        else {
                                            PolyAInfo tInfo = polyAInfoMap.get(r.tName);
                                            if (tInfo == null || isTargetPolyATContained(r, tInfo)) {
                                                containedSet.add(r.tName);
                                            }
                                        }
                                        break;
                                    case BOTH:
                                        if (isDoveTail) {
                                            containedDoveTailOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        else {
                                            if (r.qLen > r.tLen) {
                                                PolyAInfo tInfo = polyAInfoMap.get(r.tName);
                                                if (tInfo == null || isTargetPolyATContained(r, tInfo)) {
                                                    containedSet.add(r.tName);
                                                }
                                            }
                                            else if (r.tLen > r.qLen) {
                                                PolyAInfo qInfo = polyAInfoMap.get(r.qName);
                                                if (qInfo == null || isQueryPolyATContained(r, qInfo)) {
                                                    containedSet.add(r.qName);
                                                    prevContained = true;
                                                }
                                            }
                                            else {
                                                if (r.qName.compareTo(r.tName) > 0) {
                                                    PolyAInfo tInfo = polyAInfoMap.get(r.tName);
                                                    if (tInfo == null || isTargetPolyATContained(r, tInfo)) {
                                                        containedSet.add(r.tName);
                                                    }
                                                }
                                                else {
                                                    PolyAInfo qInfo = polyAInfoMap.get(r.qName);
                                                    if (qInfo == null || isQueryPolyATContained(r, qInfo)) {
                                                        containedSet.add(r.qName);
                                                        prevContained = true;
                                                    }
                                                }
                                            }
                                        }
                                        break;
                                    case NEITHER:
                                        if (isDoveTail) {
                                            pendingOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        else {
                                            nonStandardOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        break;
                                }
                            }
                        }
                    }
                }
            }
            else {
                for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
                    reader.next(r);
                    if (!r.qName.equals(r.tName)) {
                        if ((hasLargeOverlap(r) || isContainmentPafRecord(r)) && hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {                        
                            if (prevName == null) {
                                prevName = r.qName;
                            }
                            else if (!prevName.equals(r.qName)) {
                                if (!pendingOverlaps.isEmpty()) { 
                                    if (!prevContained) {
                                        overlaps.addAll(pendingOverlaps);
                                    }
                                    pendingOverlaps.clear();
                                }

                                if (this.cutRevCompArtifact && !nonStandardOverlaps.isEmpty()) {
                                    if (!prevContained) {
                                        // remove artifacts based on non standard overlaps
                                        nonStandardOverlaps.sort(overlapComparator);
                                        removeReverseComplementArtifacts(nonStandardOverlaps, containedSet);
                                    }
                                    nonStandardOverlaps.clear();
                                }

                                prevName = r.qName;
                                prevContained = containedSet.contains(prevName);
                            }

                            if (!prevContained && !containedSet.contains(r.tName)) {
                                CONTAIN_STATUS c = getContained(r);
                                boolean isDoveTail = isDovetailPafRecord(r);
                                switch (c) {
                                    case QUERY:
                                        if (isDoveTail) {
                                            containedDoveTailOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        else {
                                            containedSet.add(r.qName);
                                            prevContained = true;
                                        }
                                        break;
                                    case TARGET:
                                        if (isDoveTail) {
                                            containedDoveTailOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        else {
                                            containedSet.add(r.tName);
                                        }
                                        break;
                                    case BOTH:
                                        if (isDoveTail) {
                                            containedDoveTailOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        else {
                                            if (r.qLen > r.tLen) {
                                                containedSet.add(r.tName);
                                            }
                                            else if (r.tLen > r.qLen) {
                                                containedSet.add(r.qName);
                                                prevContained = true;
                                            }
                                            else {
                                                if (r.qName.compareTo(r.tName) > 0) {
                                                    containedSet.add(r.tName);
                                                }
                                                else {
                                                    containedSet.add(r.qName);
                                                    prevContained = true;
                                                }
                                            }
                                        }
                                        break;
                                    case NEITHER:
                                        if (isDoveTail) {
                                            pendingOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        else {
                                            nonStandardOverlaps.add(pafToStrandedOverlap(r));
                                        }
                                        break;
                                }
                            }
                        }
                    }
                }
            }
            
            reader.close();

            if (prevName != null && !pendingOverlaps.isEmpty()) {
                if (!prevContained) {
                    overlaps.addAll(pendingOverlaps);
                }
                pendingOverlaps.clear();
            }
            
            if (this.cutRevCompArtifact && !nonStandardOverlaps.isEmpty()) {
                if (!prevContained) {
                    // remove artifacts based on non standard overlaps
                    nonStandardOverlaps.sort(overlapComparator);
                    removeReverseComplementArtifacts(nonStandardOverlaps, containedSet);
                }
                nonStandardOverlaps.clear();
            }
                        
            if (polyAInfoMap != null) {
                for (StrandedOverlap r : overlaps) {
                    if (!containedSet.contains(r.qName) &&
                        !containedSet.contains(r.tName)) {
                        PolyAInfo qInfo = polyAInfoMap.get(r.qName);
                        PolyAInfo tInfo = polyAInfoMap.get(r.tName);
                        if (qInfo != null || tInfo != null) {
                            addEdges(r, qInfo, tInfo);
                        }
                        else {
                            addEdges(r);
                        }
                    }
                }
                
                for (StrandedOverlap r : containedDoveTailOverlaps) {
                    String tName = r.tName;
                    String qName = r.qName;
                    if (!containedSet.contains(qName) &&
                        !containedSet.contains(tName)) {
                        PolyAInfo qInfo = polyAInfoMap.get(r.qName);
                        PolyAInfo tInfo = polyAInfoMap.get(r.tName);
                        if (qInfo != null || tInfo != null) {
                            addEdges(r, qInfo, tInfo);
                        }
                        else {
                            addEdges(r);
                        }
                    }
                }

                // remove reverse-complement vertices of polyA tails and vertices containing polyT heads
                for (Map.Entry<String, PolyAInfo> e : polyAInfoMap.entrySet()) {
                    String name = e.getKey();
                    PolyAInfo info = e.getValue();
                    Interval polyAInterval = info.polyATail;
                    Interval polyTInterval = info.polyTHead;

                    if (polyAInterval != null && polyTInterval == null) {
                        String vid = name + '+';
                        if (graph.containsVertex(vid) && graph.outDegreeOf(vid) == 0) {
                            ArrayDeque<String> toBeRemoved = new ArrayDeque<>();
                            toBeRemoved.add(name + '-');

                            // look for reads overlapping this polyA tail
                            int polyAStart = polyAInterval.start;
                            int polyAEnd = polyAInterval.end;
                            for (OverlapEdge edge : graph.incomingEdgesOf(vid)) {
                                if (getOverlap(polyAStart, polyAEnd, edge.sinkStart, edge.sinkEnd) > 0) {
                                    toBeRemoved.add(getReverseComplementID(graph.getEdgeSource(edge)));
                                }
                            }

                            graph.removeAllVertices(toBeRemoved);
                        }
                    }
                    else if (polyAInterval == null && polyTInterval != null) {
                        String vid = name + '+';
                        if (graph.containsVertex(vid) && graph.inDegreeOf(vid) == 0) {
                            ArrayDeque<String> toBeRemoved = new ArrayDeque<>();
                            toBeRemoved.add(vid);

                            // look for reads overlapping this polyT head
                            int polyTStart = polyTInterval.start;
                            int polyTEnd = polyTInterval.end;
                            for (OverlapEdge edge : graph.outgoingEdgesOf(vid)) {
                                if (getOverlap(polyTStart, polyTEnd, edge.sourceStart, edge.sourceEnd) > 0) {
                                    toBeRemoved.add(graph.getEdgeTarget(edge));
                                }
                            }

                            graph.removeAllVertices(toBeRemoved);
                        }
                    }
                }
            }
            else {
                for (StrandedOverlap r : overlaps) {
                    if (!containedSet.contains(r.qName) &&
                        !containedSet.contains(r.tName)) {
                        addEdges(r);
                    }
                }
                
                for (StrandedOverlap r : containedDoveTailOverlaps) {
                    String tName = r.tName;
                    String qName = r.qName;
                    if (!containedSet.contains(qName) &&
                        !containedSet.contains(tName)) {
                        addEdges(r);
                    }
                }
            }
        }
        
        printMessage("Parsed " + NumberFormat.getInstance().format(reader.getNumRecords()) + " overlap records in " + timer.elapsedDHMS());
        
        return containedSet;
    }

    private boolean isQueryPolyATContained(OverlapCoords r, PolyAInfo qInfo) {
        return (qInfo.polyATail == null || 
                r.qEnd >= Math.max(qInfo.polyATail.start + minPolyALength, qInfo.polyATail.end - minPolyALength)) &&
                (qInfo.polyTHead == null ||
                r.qStart <= Math.min(qInfo.polyTHead.start + minPolyALength, qInfo.polyTHead.end - minPolyALength));
    }
    
    private boolean isTargetPolyATContained(OverlapCoords r, PolyAInfo tInfo) {
        return (tInfo.polyATail == null || 
                r.tEnd >= Math.max(tInfo.polyATail.end - minPolyALength, tInfo.polyATail.end - minPolyALength)) &&
                (tInfo.polyTHead == null || 
                r.tStart <= Math.min(tInfo.polyTHead.start + minPolyALength, tInfo.polyTHead.end - minPolyALength));
    }
    
    public void extractSimplePaths(String outFastaPath) throws IOException {
        HashSet<String> containedSet = populateGraphFromOverlaps();
        
        Timer timer = new Timer();
        
        if (!containedSet.isEmpty()) {
            printMessage("contained reads: " + NumberFormat.getInstance().format(containedSet.size()));
        }
        
        Set<String> vertexSet = graph.vertexSet();
        if (!vertexSet.isEmpty()) {
            int numDovetailReads = vertexSet.size();
            if (!stranded) {
                numDovetailReads /= 2;
            }
            printMessage("dovetail reads:  " + NumberFormat.getInstance().format(numDovetailReads));
        }
        
        Set<OverlapEdge> edgeSet = graph.edgeSet();
        
        if (edgeSet.size() > 1) {
            //printGraph();
            printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));
            
            ArrayDeque<String> removedVertexes = removeRedundantNodes();
            if (stranded) {
                for (String vid : removedVertexes) {
                    containedSet.add(getVertexName(vid));
                }
            }
            else {
                for (String vid : removedVertexes) {
                    if (!graph.containsVertex(getReverseComplementID(vid))) {
                        containedSet.add(getVertexName(vid));
                    }
                }
            }
            
            removeTransitiveEdges();
            
            //printGraph();
            printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));
            writeGraph(this.seqFastaPath + ".dot.gz");
        }
        
        HashMap<String, BitSequence> dovetailReadSeqs = new HashMap<>(vertexSet.size());
        FastaReader fr = new FastaReader(seqFastaPath);
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        long seqID = 0;
        long originalNumSeq = 0;
        while (fr.hasNext()) {
            ++originalNumSeq;
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            if (vertexSet.contains(name + '+') || (!stranded && vertexSet.contains(name + '-'))) {
                dovetailReadSeqs.put(name, new BitSequence(nameSeq[1]));
            }
            else if (!containedSet.contains(name)) {
                // this sequence either contains shorter sequences or it has no overlaps with others
                fw.write("unitig" + Long.toString(++seqID), nameSeq[1]);
            }
        }
        fr.close();
        
        // extract paths
        HashSet<String> visited = new HashSet<>();
        
        if (!stranded) {
            vertexSet = graph.vertexSet();
            // prioritize extension of vertexes with known strand
            for (String seedVid : vertexSet) {
                if (!vertexSet.contains(getReverseComplementID(seedVid))) {
                    String name = getVertexName(seedVid);
                    if (!visited.contains(name)) {
                        ArrayDeque<String> path = getUnambiguousExtension(seedVid);

                        String header = "unitig" + Long.toString(++seqID);
                        if (path.size() > 1) {
                            header += " path=[" + String.join(",", path) + "]";
                        }

                        fw.write(header, assemblePath(path, dovetailReadSeqs));

                        for (String vid : path) {
                            visited.add(getVertexName(vid));
                        }
                    }
                }
            }
        }
        
        for (String seedVid : graph.vertexSet()) {
            String name = getVertexName(seedVid);
            if (!visited.contains(name)) {
                ArrayDeque<String> path = getUnambiguousExtension(seedVid);
                
                String header = "unitig" + Long.toString(++seqID);
                if (path.size() > 1) {
                    header += " path=[" + String.join(",", path) + "]";
                }
                
                fw.write(header, assemblePath(path, dovetailReadSeqs));

                for (String vid : path) {
                    visited.add(getVertexName(vid));
                }
            }
        }
        
        fw.close();
        
        printMessage("before: " + NumberFormat.getInstance().format(originalNumSeq) + "\tafter: " + NumberFormat.getInstance().format(seqID));
        printMessage("Laid out paths in " + timer.elapsedDHMS());
    }
    
    public int removeDifferentMagnitudeEdges(HashMap<String, Float> readCounts, float maxLogDiff) {
        int numEdgesRemoved = 0;
        for (Entry<String, Float> entry : readCounts.entrySet()) {
            String name = entry.getKey();
            float count = entry.getValue();
            String vid = name + '+';
            String vidRC = name + '-';
            ArrayDeque<OverlapEdge> toBeRemoved = new ArrayDeque();
            
            if (graph.containsVertex(vid)) {
                for (OverlapEdge e : graph.incomingEdgesOf(vid)) {
                    String pName = getVertexName(graph.getEdgeSource(e));
                    Float pCount = readCounts.get(pName);
                    if (pCount == null ||
                            (count > pCount && Math.log10(pCount) + maxLogDiff < Math.log10(count)) ||
                            (count <= pCount && Math.log10(count) + maxLogDiff < Math.log10(pCount))) {
                        toBeRemoved.add(e);
                    }
                }
                
                for (OverlapEdge e : graph.outgoingEdgesOf(vid)) {
                    String pName = getVertexName(graph.getEdgeTarget(e));
                    Float pCount = readCounts.get(pName);
                    if (pCount == null ||
                            (count > pCount && Math.log10(pCount) + maxLogDiff < Math.log10(count)) ||
                            (count <= pCount && Math.log10(count) + maxLogDiff < Math.log10(pCount))) {
                        toBeRemoved.add(e);
                    }
                }
            }
            
            if (graph.containsVertex(vidRC)) {
                for (OverlapEdge e : graph.incomingEdgesOf(vidRC)) {
                    String pName = getVertexName(graph.getEdgeSource(e));
                    Float pCount = readCounts.get(pName);
                    if (pCount == null ||
                            (count > pCount && Math.log10(pCount) + maxLogDiff < Math.log10(count)) ||
                            (count <= pCount && Math.log10(count) + maxLogDiff < Math.log10(pCount))) {
                        toBeRemoved.add(e);
                    }
                }
                
                for (OverlapEdge e : graph.outgoingEdgesOf(vidRC)) {
                    String pName = getVertexName(graph.getEdgeTarget(e));
                    Float pCount = readCounts.get(pName);
                    if (pCount == null ||
                            (count > pCount && Math.log10(pCount) + maxLogDiff < Math.log10(count)) ||
                            (count <= pCount && Math.log10(count) + maxLogDiff < Math.log10(pCount))) {
                        toBeRemoved.add(e);
                    }
                }
            }
            
            numEdgesRemoved += toBeRemoved.size();
            graph.removeAllEdges(toBeRemoved);
        }
        
        return numEdgesRemoved;
    }

    private static class PolyAScore {
        public float forwardScore = 0;
        public float reverseScore = 0;
    } 
    
    private HashMap<String, PolyAScore> getPolyAScores(String pafPath,
            String polyAReadNamesPath,
            Set<String> skipSet) throws FileNotFoundException, IOException {
        
        // store the names of vertexes containing polyA tails and are at/near a dead end of the graph
        HashSet<String> polyAVertexNames = new HashSet<>();
        Set<String> vertexSet = graph.vertexSet();
        for (String vid : vertexSet) {
            if (!graph.containsVertex(getReverseComplementID(vid))) {
                polyAVertexNames.add(getVertexName(vid));
            }
        }
        
        // store the names of reads that have a *potential* polyA tail
        HashSet<String> polyAReadNames = new HashSet<>();
        fileToStringCollection(polyAReadNamesPath, polyAReadNames);
        
        HashMap<String, PolyAScore> scores = new HashMap<>();
        
        ArrayList<String> currentTNames = new ArrayList<>();
        ArrayList<PolyAScore> currentScores = new ArrayList<>();
        PafReader reader = new PafReader(pafPath);
        String prevName = null;
        boolean hasPolyA = false;
        boolean isPolyAContained = false;
        for (PafRecord r = new PafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if (!r.qName.equals(prevName)) {
                if (hasPolyA && !currentTNames.isEmpty()) {
                    if (!isPolyAContained) {
                        int numTargets = currentTNames.size();
                        for (int i=0; i<numTargets; ++i) {
                            String tName = currentTNames.get(i);
                            PolyAScore score = currentScores.get(i);

                            PolyAScore bestScore = scores.get(tName);
                            if (bestScore == null) {
                                scores.put(tName, score);
                            }
                            else {
                                bestScore.reverseScore = Math.max(bestScore.reverseScore, score.reverseScore);
                                bestScore.forwardScore = Math.max(bestScore.forwardScore, score.forwardScore);
                            }
                        }
                    }
                    
                    currentTNames = new ArrayList<>();
                    currentScores = new ArrayList<>();
                }
                
                prevName = r.qName;
                hasPolyA = polyAReadNames.contains(r.qName);
                isPolyAContained = false;
            }
            
            if (hasPolyA && !skipSet.contains(r.tName)) {
                // the polyA tail cannot be contained in another sequence
                // except if that is also a true polyA tail
                if (isContainmentPafRecord(r) && r.qEnd >= r.qLen) {
                    // the query's polyA tail is contained
                    
                    if (polyAVertexNames.contains(r.tName)) {
                        // target has a true polyA tail
                        
                        PolyAInfo info = polyAInfoMap.get(r.tName);
                        if (info != null && info.polyATail != null) {
                            if (r.tEnd < info.polyATail.start) {
                                // the query's polyA aligns upstream of the target's polyA tail
                                // therefore, query is not a true polyA tail
                                isPolyAContained = true;
                            }
                        }
                    }
                    else {
                        if (vertexSet.contains(r.tName + '+')) {
                            // target is in the graph, but it does not have a polyA tail
                            isPolyAContained = true;
                        }
                        else {
                            // target is an orphan
                            PolyAInfo info = polyAInfoMap.get(r.tName);
                            if (info == null || /* target has no polyA info */
                                    info.polyATail == null || /* target has no polyA tail */
                                    info.polyTHead != null /* target has both polyA tail and polyT head */) {
                                isPolyAContained = true;
                            }
                        }
                    }
                }

                if (isQueryEdgeSink(r, maxEdgeClip)) {
                    PolyAScore score = new PolyAScore();
                    
                    float alignedProportion = (r.tEnd - r.tStart)/(float)r.tLen;
                    
                    if (r.reverseComplemented) {
                        score.reverseScore = alignedProportion;
                    }
                    else {
                        score.forwardScore = alignedProportion;
                    }
                }
            }
        }
        reader.close();
        
        return scores;
    }
    
    private void pruneGraphWithPolyAInfo(String mappingPafPath, String polyAReadNamesPath, Set<String> skipSet) throws IOException {        
        HashMap<String, PolyAScore> polyAScores = getPolyAScores(mappingPafPath, polyAReadNamesPath, skipSet);
        
        for (Entry<String, PolyAScore> e : polyAScores.entrySet()) {
            String name = e.getKey();
            PolyAScore score = e.getValue();
            if (score.forwardScore > 0 && score.reverseScore == 0) {
                // remove minus vertex
                if (graph.containsVertex(name + '+')) {
                    graph.removeVertex(name + '-');
                }
            }
            else if (score.forwardScore == 0 && score.reverseScore > 0) {
                // remove plus vertex
                if (graph.containsVertex(name + '-')) {
                    graph.removeVertex(name + '+');
                }
            }
        }

        ArrayDeque<OverlapEdge> toBeRemoved = new ArrayDeque<>();
        for (String vid : graph.vertexSet()) {
            if (!graph.containsVertex(getReverseComplementID(vid))) {
                for (OverlapEdge edge : graph.outgoingEdgesOf(vid)) {
                    String sid = graph.getEdgeTarget(edge);
                    if (graph.containsVertex(getReverseComplementID(sid))) {
                        toBeRemoved.add(edge);
                    }
                }
            }
        }
        graph.removeAllEdges(toBeRemoved);
    }
    
    private void filterEdges(HashMap<String, Float> readCounts, String sampleReadLengthsPath) throws IOException {
        int[] sampleLengths = readIntArrayFromFile(sampleReadLengthsPath);
        int maxLength = sampleLengths[sampleLengths.length-1];
        EmpiricalDistribution readLenDist = EmpiricalDistribution.fit(sampleLengths, IntSet.of(maxLength+1));
        
        ArrayDeque<OverlapEdge> toBeRemoved = new ArrayDeque<>();
        int numSupported = 0;
        for (OverlapEdge e : graph.edgeSet()) {
            float overlapSize = ((e.sinkEnd - e.sinkStart) + (e.sourceEnd - e.sourceStart))/2f;
            double supportingReads = graph.getEdgeWeight(e) - DEFAULT_EDGE_WEIGHT;
            if (supportingReads > 0) {
                ++numSupported;
            } 
            
            if (overlapSize < maxLength) {
                String sID = graph.getEdgeSource(e);
                String tID = graph.getEdgeTarget(e);

                String sName = getVertexName(sID);
                String tName = getVertexName(tID);

                Float sCount = readCounts.get(sName);
                Float tCount = readCounts.get(tName);

                if (sCount == null ) {
                    sCount = (float) 0;
                }

                if (tCount == null ) {
                    tCount = (float) 0;
                }

                // probability of finding one read smaller than the overlap
                double p = readLenDist.cdf(overlapSize);

                // expected number of reads spanning the overlap
                double c = Math.floor(Math.max(sCount, tCount));
                
                if (supportingReads < c) {
                    // binomial test
                    BinomialDistribution bd = new BinomialDistribution((int)c, 1-p);
                    if (bd.cdf(supportingReads) < 0.001) {
                        toBeRemoved.add(e);
                    }
                }
            }
        }
        
        System.out.println("Supported edges: " + numSupported);
        
        graph.removeAllEdges(toBeRemoved);
    }
    
    public void extractGreedyPaths(String outFastaPath, String mappingPafPath,
            String polyAReadNamesPath, String sampleReadLengthsPath,
            String txptNamePrefix, boolean writeUracil) throws IOException {
        HashSet<String> containedSet = populateGraphFromOverlaps();
        
        if (!containedSet.isEmpty()) {
            printMessage("contained reads: " + NumberFormat.getInstance().format(containedSet.size()));
        }
        
        Set<String> vertexSet = graph.vertexSet();
        if (!vertexSet.isEmpty()) {
            int numDovetailReads = vertexSet.size();
            if (!stranded) {
                numDovetailReads /= 2;
            }
            printMessage("dovetail reads:  " + NumberFormat.getInstance().format(numDovetailReads));
        }
        
        Set<OverlapEdge> edgeSet = graph.edgeSet();
        
        if (edgeSet.size() > 1) {            
            printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));
            printMessage("Removing redundant vertexes...");
            ArrayDeque<String> removedVertexes = removeRedundantNodes();
            if (stranded) {
                for (String vid : removedVertexes) {
                    containedSet.add(getVertexName(vid));
                }
            }
            else {
                for (String vid : removedVertexes) {
                    if (!graph.containsVertex(getReverseComplementID(vid))) {
                        containedSet.add(getVertexName(vid));
                    }
                }
            }
            
            printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));
            printMessage("Removing transitive edges...");
            removeTransitiveEdges();
            
            if (!stranded && polyAReadNamesPath != null) {
                printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));
                printMessage("Pruning graph with polyA information...");
                pruneGraphWithPolyAInfo(mappingPafPath, polyAReadNamesPath, containedSet);
            }       
            
            //printGraph();
            printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));
            writeGraph(seqFastaPath + ".dot.gz");
        }
                
        printMessage("Tallying read counts...");

        Timer timer = new Timer();
        //HashMap<String, Float> readCounts = getReadCounts(mappingPafPath, vertexNames, maxEdgeClip);
        HashMap<String, Float> readCounts = getLengthNormalizedReadCounts(mappingPafPath, containedSet);
        printMessage("Counts tallied for " + readCounts.size() + " sequences in " + timer.elapsedDHMS());
        
        //printCounts(readCounts);
        
        printMessage("Pruning graph with read count information...");
        filterEdges(readCounts, sampleReadLengthsPath);
        printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));
        
//        // assign read count to edge
//        printMessage("Adding read counts to graph...");
//        timer.start();
//        edgeSet = graph.edgeSet();
//        for (OverlapEdge e : edgeSet) {
//            String source = getVertexName(graph.getEdgeSource(e));
//            String target = getVertexName(graph.getEdgeTarget(e));
////            Float sCount = readCounts.get(source.substring(0, source.length()-1));
////            Float tCount = readCounts.get(target.substring(0, target.length()-1));
//            Float sCount = readCounts.get(source);
//            Float tCount = readCounts.get(target);
//            float w = 0;
//            if (sCount != null && tCount != null) {
//                w = Math.min(sCount, tCount);
//            }
//            else if (sCount != null) {
//                w = sCount;
//            }
//            else if (tCount != null) {
//                w = tCount;
//            }
//            graph.setEdgeWeight(e, w);
//        }
//        
//        printMessage("Counts added in " + timer.elapsedDHMS());
        printMessage("Extracting vertex sequences...");

        timer.start();
        HashMap<String, BitSequence> dovetailReadSeqs = new HashMap<>(vertexSet.size());
        FastaReader fr = new FastaReader(seqFastaPath);
        FastaWriter fw = new FastaWriter(outFastaPath, false, writeUracil);
        long seqID = 0;
        long originalNumSeq = 0;
        while (fr.hasNext()) {
            ++originalNumSeq;
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            if (vertexSet.contains(name + '+') || vertexSet.contains(name + '-')) {
                dovetailReadSeqs.put(name, new BitSequence(nameSeq[1]));
            }
            else if (!containedSet.contains(name)) {
                // this sequence either contains shorter sequences or it has no overlaps with others
//                fw.write(txptNamePrefix + Long.toString(++seqID) + " l=" + nameSeq[1].length() + " c=" + readCounts.get(name) + " s=" + name, nameSeq[1]);
//                readCounts.remove(name);

                Float c = readCounts.get(name);
                if (c != null) {
                    fw.write(txptNamePrefix + Long.toString(++seqID) + " l=" + nameSeq[1].length() + " c=" + c.toString() + " s=" + name, nameSeq[1]);
                    readCounts.remove(name);
                }
            }
        }
        fr.close();
        
        printMessage("Sequences extracted in " + timer.elapsedDHMS());
        printMessage("Extracting paths...");
        
        // extract paths
        timer.start();
        HashSet<String> visited = new HashSet<>();
        for (Iterator<Entry<String, Float>> itr = readCounts.entrySet().stream().
                sorted(Entry.<String, Float>comparingByValue().reversed()).iterator();
                itr.hasNext();) {
            Entry<String, Float> e = itr.next();
            String seed = e.getKey();// + "r";
            if (!visited.contains(seed)) {
                String seedVid = seed + '+';
                if (!graph.containsVertex(seedVid)) {
                    if (stranded) {
                        continue;
                    }
                    else {
                        seedVid = seed + '-';
                        if (!graph.containsVertex(seedVid)) {
                            continue;
                        }
                    }
                }
                
                ArrayDeque<String> path = getMaxWeightExtension(seedVid, readCounts);
                
                String header = txptNamePrefix + Long.toString(++seqID);
                String seq = assemblePath(path, dovetailReadSeqs);

                header += " l=" + seq.length() + " c=" + getMinAndDecrementWeights(path, readCounts);

                if (path.size() > 1) {
                    header += " path=[" + String.join(",", path) + "]";
                }
                else {
                    header += " s=" + seed;
                }

                fw.write(header, seq);

                for (String vid : path) {
                    visited.add(getVertexName(vid));
                }
            }
        }
        
        fw.close();
        
        printMessage("before: " + NumberFormat.getInstance().format(originalNumSeq) + "\tafter: " + NumberFormat.getInstance().format(seqID));
        printMessage("Laid out paths in " + timer.elapsedDHMS());
    }

    private String getMaxWeightPredecessor(String vid, HashMap<String, Float> weights) {
        float maxWeight = 0; //Double.MIN_VALUE;
        String maxWeightPredecessor = null;
        for (String p : Graphs.predecessorListOf(graph, vid)) {
            Float w = weights.get(getVertexName(p));
            if (w != null && w > maxWeight) {
                maxWeight = w;
                maxWeightPredecessor = p;
            }
        }
        
        return maxWeightPredecessor;
    }
    
    private String getMaxWeightSuccessor(String vid, HashMap<String, Float> weights) {
        float maxWeight = 0; //Double.MIN_VALUE;
        String maxWeightSuccessor = null;
        for (String s : Graphs.successorListOf(graph, vid)) {
            Float w = weights.get(getVertexName(s));
            if (w != null && w > maxWeight) {
                maxWeight = w;
                maxWeightSuccessor = s;
            }
        }
        
        return maxWeightSuccessor;
    }
    
    private ArrayDeque<String> getGreedyExtension(String seedVid) {
        HashSet<String> visited = new HashSet<>();
        visited.add(seedVid);
        
        ArrayDeque<String> path = new ArrayDeque<>();
                
        // extend left
        String cursor = seedVid;
        while ((cursor = getClosestPredecessor(cursor)) != null) {
            if (visited.contains(cursor)) {
                break;
            }
            path.addFirst(cursor);
            visited.add(cursor);
        }
        
        path.add(seedVid);
        
        // extend right
        cursor = seedVid;
        while ((cursor = getClosestSuccessor(cursor)) != null) {
            if (visited.contains(cursor)) {
                break;
            }
            path.add(cursor);
            visited.add(cursor);
        }
        
        return path;
    }
    
    private String getClosestPredecessor(String vid) {
        int maxOverlapSize = Integer.MIN_VALUE;
        OverlapEdge closestEdge = null;
        for (OverlapEdge e : graph.incomingEdgesOf(vid)) {
            int o = getOverlapSize(e);
            if (o > maxOverlapSize) {
                maxOverlapSize = o;
                closestEdge = e;
            }
        }
        
        if (closestEdge != null) {
            return graph.getEdgeSource(closestEdge);
        }
        return null;
    }
    
    private String getClosestSuccessor(String vid) {
        int maxOverlapSize = Integer.MIN_VALUE;
        OverlapEdge closestEdge = null;
        for (OverlapEdge e : graph.outgoingEdgesOf(vid)) {
            int o = getOverlapSize(e);
            if (o > maxOverlapSize) {
                maxOverlapSize = o;
                closestEdge = e;
            }
        }
        
        if (closestEdge != null) {
            return graph.getEdgeTarget(closestEdge);
        }
        return null;
    }
    
    private int getOverlapSize(OverlapEdge e) {
        return ((e.sinkEnd - e.sinkStart) + (e.sourceEnd - e.sourceStart))/2;
    }
    
    private ArrayDeque<String> getMaxWeightExtension(String seedVid, HashMap<String, Float> weights) {
        HashSet<String> visited = new HashSet<>();
        visited.add(seedVid);
        
        ArrayDeque<String> path = new ArrayDeque<>();
        
        Float seedWeight = weights.get(getVertexName(seedVid));
        if (seedWeight != null) {
            // extend left
            String cursor = seedVid;
            while ((cursor = getMaxWeightPredecessor(cursor, weights)) != null) {
                if (visited.contains(cursor)) {
                    break;
                }
                Float weight = weights.get(getVertexName(cursor));
                if (weight != null) {
                    path.addFirst(cursor);
                    visited.add(cursor);
                }
                else {
                    break;
                }
            }

            path.add(seedVid);

            // extend right
            cursor = seedVid;
            while ((cursor = getMaxWeightSuccessor(cursor, weights)) != null) {
                if (visited.contains(cursor)) {
                    break;
                }
                Float weight = weights.get(getVertexName(cursor));
                if (weight != null) {
                    path.add(cursor);
                    visited.add(cursor);
                }
                else {
                    break;
                }
            }
        }
        
        return path;
    }
    
    private float getMinAndDecrementWeights(ArrayDeque<String> path, HashMap<String, Float> weights) {
        float minWeight = Float.MAX_VALUE;
        for (String vid : path) {
            String name = getVertexName(vid);
            float w = weights.get(name);
            if (w < minWeight) {
                minWeight = w;
            }
        }
        
        for (String vid : path) {
            String name = getVertexName(vid);
            float w = weights.get(name);
            weights.put(name, w - minWeight);
        }
        
        return minWeight;
    }
    
    public boolean layoutBackbones(String outFastaPath) throws IOException {
        HashSet<String> containedSet = populateGraphFromOverlaps();
        
        if (!containedSet.isEmpty()) {
            printMessage("contained reads: " + NumberFormat.getInstance().format(containedSet.size()));
        }
        
        Set<String> vertexSet = graph.vertexSet();
        if (!vertexSet.isEmpty()) {
            int numDovetailReads = vertexSet.size();
            if (!stranded) {
                numDovetailReads = numDovetailReads/2;
            }
            printMessage("dovetail reads:  " + NumberFormat.getInstance().format(numDovetailReads));
        }
        
        Set<OverlapEdge> edgeSet = graph.edgeSet();
        
        if (edgeSet.size() > 1) {
            printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));

//            ArrayDeque<String> redundantNodeNames = removeRedundantNodes();
//            for (String n : redundantNodeNames) {
//                containedSet.add(getVertexName(n));
//            }
//            
//            printGraph();
//            numEdges = graph.edgeSet().size();
//            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));

            resolveJunctions();
            printMessage("G: |V|=" + NumberFormat.getInstance().format(vertexSet.size()) + " |E|=" + NumberFormat.getInstance().format(edgeSet.size()));
        }
        
        // extract read sequences
        int numVertexReads = vertexSet.size();
        if (!stranded) {
            numVertexReads = numVertexReads/2;
        }
        
        HashMap<String, BitSequence> dovetailReadSeqs = new HashMap<>(numVertexReads);
        FastaReader fr = new FastaReader(seqFastaPath);
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        long seqID = 0;
        long originalNumSeq = 0;
        while (fr.hasNext()) {
            ++originalNumSeq;
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            if (vertexSet.contains(name + '+')) {
                dovetailReadSeqs.put(name, new BitSequence(nameSeq[1]));
            }
            else if (!containedSet.contains(name)) {
                // this sequence either contains shorter sequences or it has no overlaps with others
                fw.write(Long.toString(++seqID), nameSeq[1]);
            }
        }
        fr.close();
        
        // layout unambiguous paths
        HashSet<String> visitedReadNames = new HashSet<>(numVertexReads);
        for (String n : vertexSet) {
            if (!visitedReadNames.contains(getVertexName(n))) {
                
                ArrayDeque<String> path = getUnambiguousLeftExtension(n);
                path.add(n);
                if (!graph.containsEdge(n, path.getFirst())) {
                    // detect cycles
                    path.addAll(getUnambiguousRightExtension(n));
                }
                
                String backbone = assemblePath(path, dovetailReadSeqs);
                
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
        
        printMessage("before: " + NumberFormat.getInstance().format(originalNumSeq) + "\tafter: " + NumberFormat.getInstance().format(seqID));
        
        return originalNumSeq > seqID;
    }

//    private class OverlapCoords {
//        int qStart, qEnd, tStart, tEnd;
//    }
    
    private class TargetOverlap extends OverlapCoords {
        String tName;
    }
    
    private class QueryOverlap extends OverlapCoords {
        String qName;
    }
    
    private class Overlap extends OverlapCoords {
        String qName, tName;
        int qLen, tLen;
    }
    
    private class StrandedOverlap extends Overlap {
        boolean reverseComplemented;
    }
    
    private Overlap pafToOverlap(PafRecord r) {
        Overlap o = new Overlap();
        o.qName = r.qName;
        o.tName = r.tName;
        o.qLen = r.qLen;
        o.tLen = r.tLen;
        o.qStart = r.qStart;
        o.qEnd = r.qEnd;
        o.tStart = r.tStart;
        o.tEnd = r.tEnd;
        return o;
    }
    
    private StrandedOverlap pafToStrandedOverlap(PafRecord r) {
        StrandedOverlap o = new StrandedOverlap();
        o.qName = r.qName;
        o.tName = r.tName;
        o.qLen = r.qLen;
        o.tLen = r.tLen;
        o.qStart = r.qStart;
        o.qEnd = r.qEnd;
        o.tStart = r.tStart;
        o.tEnd = r.tEnd;
        o.reverseComplemented = r.reverseComplemented;
        return o;
    }
    
    private QueryOverlap pafToQueryOverlap(PafRecord r) {
        QueryOverlap o = new QueryOverlap();
        o.qName = r.qName;
        o.qStart = r.qStart;
        o.qEnd = r.qEnd;
        o.tStart = r.tStart;
        o.tEnd = r.tEnd;
        return o;
    }
    
    private TargetOverlap pafToTargetOverlap(PafRecord r) {
        TargetOverlap o = new TargetOverlap();
        o.tName = r.tName;
        o.qStart = r.qStart;
        o.qEnd = r.qEnd;
        o.tStart = r.tStart;
        o.tEnd = r.tEnd;
        return o;
    }
    
    private boolean layoutStrandedBackbones(String outFastaPath) throws IOException {
        HashSet<String> containedSet = populateGraphFromOverlaps();
        
        if (!containedSet.isEmpty()) {
            printMessage("contained reads: " + NumberFormat.getInstance().format(containedSet.size()));
        }
        
        Set<String> vertexSet = graph.vertexSet();
        if (!vertexSet.isEmpty()) {
            printMessage("dovetail reads:  " + NumberFormat.getInstance().format(vertexSet.size()));
        }
        
        int numEdges = graph.edgeSet().size();        
        if (numEdges > 1) {
            printMessage("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
            
//            ArrayDeque<String> redundantNodeNames = removeRedundantNodes();
//            for (String n : redundantNodeNames) {
//                containedSet.add(getVertexName(n));
//            }
//            numEdges = graph.edgeSet().size();
//            System.out.println("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));           
            
            resolveJunctions();
            numEdges = graph.edgeSet().size();
            printMessage("G: |V|=" + NumberFormat.getInstance().format(graph.vertexSet().size()) + " |E|=" + NumberFormat.getInstance().format(numEdges));
        }
        
        // extract read sequences
        HashMap<String, BitSequence> dovetailReadSeqs = new HashMap<>(vertexSet.size());
        FastaReader fr = new FastaReader(seqFastaPath);
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        long seqID = 0;
        long originalNumSeq = 0;
        while (fr.hasNext()) {
            ++originalNumSeq;
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            if (vertexSet.contains(name + '+')) {
                String seq = nameSeq[1];
                dovetailReadSeqs.put(name, new BitSequence(seq));
            }
            else if (!containedSet.contains(name)) {
                // this sequence either contains shorter sequences or it has no overlaps with others
                fw.write(Long.toString(++seqID), nameSeq[1]);
            }
        }
        fr.close();
        
        // layout unambiguous paths
        HashSet<String> visitedReadNames = new HashSet<>(vertexSet.size());
        for (String n : vertexSet) {
            if (!visitedReadNames.contains(getVertexName(n))) {                
                ArrayDeque<String> path = getUnambiguousLeftExtension(n);
                path.add(n);
                if (!graph.containsEdge(n, path.getFirst())) {
                    // detect cycles
                    path.addAll(getUnambiguousRightExtension(n));
                }
                
                //print(path, ',')
                String backbone = assemblePath(path, dovetailReadSeqs);
                
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
        
        printMessage("before: " + NumberFormat.getInstance().format(originalNumSeq) + "\tafter: " + NumberFormat.getInstance().format(seqID));
        
        return originalNumSeq > seqID;
    }
    
    private class PolyAInfo {
        Interval polyATail;
        Interval polyTHead;
        boolean hasPas;
        boolean hasPasRC;
        
        public PolyAInfo(Interval polyATail, Interval polyTHead,
            boolean hasPas, boolean hasPasRC) {
            this.polyATail = polyATail;
            this.polyTHead = polyTHead;
            this.hasPas = hasPas;
            this.hasPasRC = hasPasRC;
        }
    }
    
    private HashMap<String, PolyAInfo> getPolyAInfo(PolyATailFinder finder) throws IOException {
        HashMap<String, PolyAInfo> info = new HashMap<>();
        FastaReader fr = new FastaReader(seqFastaPath);
        while (fr.hasNext()) {
            String[] nameSeq = fr.nextWithName();
            String seq = nameSeq[1];
            Interval polyA = finder.findPolyATail(seq);
            Interval polyT = finder.findPolyTHead(seq);
            if (polyA != null || polyT != null) {
                boolean hasPas = polyA == null ? false : finder.hasPolyASignal(seq, polyA.start);
                boolean hasPasRC = polyT == null ? false : finder.hasPolyASignalRC(seq, polyT.end);
                info.put(nameSeq[0], new PolyAInfo(polyA, polyT, hasPas, hasPasRC));
            }
        }
        fr.close();
        return info;
    }
    
    private static Overlap getOverlapContainer(Overlap q, Collection<Overlap> list, float maxProportion) {
        int maxOverlapLen = 0;
        Overlap container = null;
        for (Overlap other : list) {
            int overlapLen = getOverlap(q.qStart, q.qEnd, other.qStart, other.qEnd);
            if (overlapLen > maxOverlapLen) {
                maxOverlapLen = overlapLen;
                container = other;
            }
        }
        
        if (maxOverlapLen >= maxProportion * (q.qEnd - q.qStart)) {
            return container;
        }
        
        return null;
    }
    
    private void updateCounts(ArrayList<StrandedOverlap> targets, HashMap<String, Float> counts) {
        if (!targets.isEmpty()) {
            int numTargets = targets.size();
            if (numTargets == 1) {
                Overlap t = targets.get(0);
                Float c = counts.get(t.tName);
                if (c == null) {
                    c = (t.tEnd - t.tStart)/(float) t.tLen;
                }
                else {
                    c += (t.tEnd - t.tStart)/(float) t.tLen;
                }
                counts.put(t.tName, c);
            }
            else {
                Collections.sort(targets, qStartComparator);
                
                for (int i=0; i<numTargets; ++i) {
                    StrandedOverlap left = targets.get(i);
                    String leftVid = left.tName;
                    String leftVidRC = left.tName;
                    if (left.reverseComplemented) {
                        leftVid += '-';
                        leftVidRC += '+';
                    }
                    else {
                        leftVid += '+';
                        leftVidRC += '-';
                    }
                    
                    for (int j=i+1; j<numTargets; ++j) {
                        StrandedOverlap right = targets.get(j);
                        if (right.qStart > left.qEnd) {
                            break;
                        }
                        
                        if (isForwardDoveTail(left.qStart, left.qEnd, right.qStart, right.qEnd)) {
                            String rightVid = right.tName;
                            String rightVidRC = right.tName;
                            if (right.reverseComplemented) {
                                rightVid += '-';
                                rightVidRC += '+';
                            }
                            else {
                                rightVid += '+';
                                rightVidRC += '-';
                            }
                            
                            OverlapEdge edge = graph.getEdge(leftVid, rightVid);
                            if (edge != null) {
                                double weight = graph.getEdgeWeight(edge);
                                graph.setEdgeWeight(edge, weight + 1);
                            }
                            
                            edge = graph.getEdge(rightVidRC, leftVidRC);
                            if (edge != null) {
                                double weight = graph.getEdgeWeight(edge);
                                graph.setEdgeWeight(edge, weight + 1);
                            }
                        }
                    }
                }
                
                Collections.sort(targets, overlapSizeComparator);

                ArrayList<Overlap> targetsKept = new ArrayList<>();
                HashMap<Overlap, ArrayList<Overlap>> multiTargets = new HashMap<>();
                for (Overlap m : targets) {
                    Overlap c = getOverlapContainer(m, targetsKept, 0.95f);
                    if (c == null) {
                        // region not contained
                        targetsKept.add(m);
                    }
                    else if (m.qEnd - m.qStart >= (c.qEnd - c.qStart) * 0.95f) {
                        // region multimaps
                        ArrayList<Overlap> multimaps = multiTargets.get(c);
                        if (multimaps == null) {
                            multimaps = new ArrayList<>();
                            multiTargets.put(c, multimaps);
                        }
                        multimaps.add(m);
                    }
                }

                for (Overlap t : targetsKept) {
                    ArrayList<Overlap> multimaps = multiTargets.get(t);
                    if (multimaps != null) {
                        // fractional assignment of multimapped region
                        float fraction = 1f/(multimaps.size() + 1);

                        Float c = counts.get(t.tName);
                        if (c == null) {
                            c = (t.tEnd - t.tStart)/(float) t.tLen * fraction;
                        }
                        else {
                            c += (t.tEnd - t.tStart)/(float) t.tLen * fraction;
                        }
                        counts.put(t.tName, c);

                        for (Overlap mm : multimaps) {
                            c = counts.get(mm.tName);
                            if (c == null) {
                                c = (mm.tEnd - mm.tStart)/(float) mm.tLen * fraction;
                            }
                            else {
                                c += (mm.tEnd - mm.tStart)/(float) mm.tLen * fraction;
                            }
                            counts.put(mm.tName, c);
                        }
                    }
                    else {
                        Float c = counts.get(t.tName);
                        if (c == null) {
                            c = (t.tEnd - t.tStart)/(float) t.tLen;
                        }
                        else {
                            c += (t.tEnd - t.tStart)/(float) t.tLen;
                        }
                        counts.put(t.tName, c);
                    }
                }
            }
        }
    }
    
    public HashMap<String, Float> getLengthNormalizedReadCounts(String pafPath, Set<String> skipSet) throws IOException {
        HashMap<String, Float> counts = new HashMap<>();
                
        PafReader reader = new PafReader(pafPath);
        String prevName = null;
        ArrayList<StrandedOverlap> targets = new ArrayList<>();
        for (PafRecord r = new PafRecord(); reader.hasNext();) {
            reader.next(r);
            
            if (!r.qName.equals(prevName)) {
                updateCounts(targets, counts);
                targets = new ArrayList<>();
                prevName = r.qName;
            }
            
            if (!skipSet.contains(r.tName)) {
                targets.add(pafToStrandedOverlap(r));
            }
        }
        reader.close();
        
        // process the last batch of records
        updateCounts(targets, counts);
        
        return counts;
    }
    
    public static void main(String[] args) {
        try {
            //debug
            String dir = "/home/gengar/test_assemblies/ont_work/complex3";
            String paf = dir + "/pol.paf";
            String in = dir + "/rnabloom.longreads.assembly4.pol.fa";
            String out = dir + "/out.fa";
            PolyATailFinder polyaFinder = new PolyATailFinder();
            polyaFinder.setProfile(PolyATailFinder.Profile.ONT);
            polyaFinder.setWindow(0);
            
            Layout layout = new Layout(in, new BufferedInputStream(new FileInputStream(paf)), false,
                    50, 0.7f, 150, 50,
                    true, 1, true, polyaFinder);
            layout.extractSimplePaths(out);
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
    }
}
