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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.text.NumberFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.TransitiveReduction;
import org.jgrapht.alg.connectivity.KosarajuStrongConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import rnabloom.util.MiniFloat;
import rnabloom.io.CompressedFastaRecord;
import static rnabloom.io.Constants.FASTA_EXT;
import rnabloom.io.ExtendedPafRecord;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaRecord;
import rnabloom.io.FastaWriter;
import rnabloom.io.PafReader;
import rnabloom.io.PafRecord;
import static rnabloom.util.Common.convertToRoundedPercent;
import static rnabloom.util.PafUtils.*;
import static rnabloom.util.SeqUtils.stringToBytes;
import static rnabloom.util.SeqUtils.bytesToString;
import static rnabloom.util.SeqUtils.reverseComplement;
import rnabloom.util.Timer;


/**
 *
 * @author Ka Ming Nip
 */
public class Layout {
    
    private DefaultDirectedGraph<String, OverlapEdge> graph;
    private InputStream overlapPafInputStream;
    private String seqFastaPath;
    private boolean stranded;
    private int maxEdgeClip = 100;
    private float minAlnId = 0.50f;
    private int maxIndelSize = 20;
    private int minOverlapMatches = 200;
    private boolean cutRevCompArtifact = false;
    private int minNumAltReads = 0;
    
    public Layout(String seqFile, InputStream overlapPafInputStream, boolean stranded, int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact, int minSeqDepth) {
        this.graph = new DefaultDirectedGraph<>(OverlapEdge.class);
        this.overlapPafInputStream = overlapPafInputStream;
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
    
    private boolean isStrandedContainmentPafRecord(PafRecord r) {
        return !r.reverseComplemented && isContainmentPafRecord(r);
    }
        
    private boolean isDovetailPafRecord(PafRecord r) {
        return rnabloom.util.PafUtils.isDovetailPafRecord(r, maxEdgeClip);
    }
    
    private boolean isStrandedDovetailPafRecord(PafRecord r) { 
        return rnabloom.util.PafUtils.isStrandedDovetailPafRecord(r, maxEdgeClip);
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
    
    /**
     * 
     * @param outFastaPath
     * @return boolean whether number of sequences have reduced
     * @throws IOException 
     */
    public boolean writeBackboneSequences(String outFastaPath) throws IOException {
        if (stranded) {
            return layoutStrandedBackbones(outFastaPath);
        }
        else {
            return layoutBackbones(outFastaPath);
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
        if (length <= 250) {
            return Math.min(25, minOverlapMatches);
        }
        else if (length <= 500) {
            return Math.min(50, minOverlapMatches);
        }
        else if (length <= 1000) {
            return Math.min(100, minOverlapMatches);
        }
        else {
            return Math.min(200, minOverlapMatches);
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
    
    private static void updateHistogram(Histogram h, int start, int end, int binSize) {
        h.minStart = Math.min(h.minStart, start);
        h.maxEnd = Math.max(h.maxEnd, end);

        if (h.bars != null) {
            if (start > 0) {
                start = (int) Math.ceil((float) start/binSize);
            }
            else {
                start = 0;
            }

            if (end < h.length) {
                end = (int) Math.floor((float) end/binSize);
            }
            else {
                end = h.bars.length;
            }

            updateHistogram(h.bars, start, end);
        }
    }
    
    private boolean isFullyCovered(Histogram h) {
        return h.minStart <= maxEdgeClip && h.maxEnd >= h.length - maxEdgeClip;
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
        short c;
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

            if (first.start < hist.minStart + binSize) {
                first.start = hist.minStart;
            }

            if (last.end > hist.length) {
                last.end = hist.length;
            }
            else if (last.end > hist.maxEnd - binSize) {
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
    
    private static enum CONTAIN_STATUS {QUERY, TARGET, NEITHER};
    
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
    
    public long extractUnique(String outFastaPath) throws IOException {
        Timer timer = new Timer();
        
        //HashMap<String, String> longestAlts = new HashMap<>(); // read id : longest read id
        Set<String> containedSet = Collections.synchronizedSet(new HashSet<>());
        HashMap<String, Histogram> histogramMap = new HashMap<>(100000);
        final int minSegmentLength = minOverlapMatches;
        
        PafReader reader = new PafReader(overlapPafInputStream);
        final boolean checkNumAltReads = minNumAltReads > 0;
        long numRecords = 0;
        
        ArrayDeque<TargetOverlap> currRecords = new ArrayDeque<>();
        String currName = null;
        Histogram currHist = null;
        int currHistBinSize = -1;
        boolean currContained = false;
        
        for (PafRecord r = new PafRecord(); reader.hasNext();) {
            ++numRecords;
            reader.next(r);
            
            if ((!stranded || !r.reverseComplemented)&&
                    !r.qName.equals(r.tName) &&
                    hasLargeOverlap(r) && hasGoodOverlap(r)) {
                
                if (!r.qName.equals(currName)) {
                    if (currName != null && currHist != null) {
                        if (!currRecords.isEmpty()) {
                            findContainedTargetOverlaps(currName, currHist, currRecords, containedSet, histogramMap);
                            currRecords.clear();
                        }

                        if (currHist.pendingQueries != null) {
                            findContainedQueryOverlaps(currName, currHist, containedSet, histogramMap);
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
                
                updateHistogram(currHist, r.qStart, r.qEnd, currHistBinSize);

                Histogram tHist = histogramMap.get(r.tName);
                int tHistBinSize = getHistogramBinSize(r.tLen);
                if (tHist == null) {
                    tHist = new Histogram(r.tLen, r.tStart, r.tEnd, tHistBinSize);
                    histogramMap.put(r.tName, tHist);
                }
                
                updateHistogram(tHist, r.tStart, r.tEnd, tHistBinSize);
                
                if (!currContained || !containedSet.contains(r.tName)) {
                    // look for containment only if the other is not already "contained"
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
            }
        }
        reader.close();
        
        if (!currRecords.isEmpty() && currHist != null) {
            findContainedTargetOverlaps(currName, currHist, currRecords, containedSet, histogramMap);
            currRecords.clear();
        }
        
        for (HashMap.Entry<String, Histogram> e : histogramMap.entrySet()) {
            Histogram h = e.getValue();
            if (h.pendingQueries != null) {
                findContainedQueryOverlaps(e.getKey(), h, containedSet, histogramMap);
                h.pendingQueries = null;
            }
        }
        
        System.out.println("Parsed " + NumberFormat.getInstance().format(numRecords) + " overlap records in " + timer.elapsedDHMS());
        
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
                                        String header = namePrefix + segmentID + " " + seqLen + ":" + i.start + "-" + i.end + " " + seqName;
                                        fw.write(header, record.seq.substring(i.start, i.end));
                                        ++segmentID;
                                    }
                                }
                                else {
                                    Interval i = spans.peekFirst();
                                    String header = Long.toString(seqID) + " " + seqLen + ":" + i.start + "-" + i.end  + " " + seqName;
                                    fw.write(header, record.seq.substring(i.start, i.end));
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
        
        System.out.println("total reads:    " + NumberFormat.getInstance().format(originalNumSeq));
        System.out.println(" - unique:      " + NumberFormat.getInstance().format(seqID) + "\t(" + convertToRoundedPercent(seqID/(float)originalNumSeq) + " %)");
        System.out.println("   - multi-seg: " + NumberFormat.getInstance().format(numMultiSeqs));
        System.out.println("Unique reads extracted in " + timer.elapsedDHMS());
        
        return seqID;
    }
    
    public int[] extractClusters(String outdir, int maxMergedClusterSize) throws IOException {
        Timer timer = new Timer();
        PafReader reader = new PafReader(overlapPafInputStream);
        
        final boolean checkNumAltReads = minNumAltReads > 0;
        
        HashMap<String, short[]> readHistogramMap = new HashMap<>();
        final int histBinSize = 25;
        final int maxBestNeighbors = 2;
        BestNeighborPairs bestNeighbors = new BestNeighborPairs(maxBestNeighbors);
//        BestNeighbor bestNeighbors = new BestNeighbor();
                
        long records = 0;
        for (PafRecord r = new PafRecord(); reader.hasNext();) {
            ++records;
//            if (++records % 1000000 == 0) {
//                System.out.println("Parsed " + NumberFormat.getInstance().format(records) + " overlap records...");
//            }
            
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
        System.out.println("Parsed " + NumberFormat.getInstance().format(records) + " overlap records in " + timer.elapsedDHMS());
        
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

            System.out.println("Effective intervals identified in " + timer.elapsedDHMS());
            System.out.println("\t- multi-segs:\t" + multiSegmentSeqs.size());
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
        System.out.println("Overlaps sorted in " + timer.elapsedDHMS());

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
        System.out.println("Clusters formed in " + timer.elapsedDHMS());
        System.out.println("\t- clusters:\t" + numClusters);

        // assign cluster IDs
        timer.start();
        clusters.assignIDs();
        int largestClusterID = clusters.getLargestClusterID();
        int largestClusterSize= clusters.getLargestClusterSize();
        System.out.println("Cluster IDs assigned in " + timer.elapsedDHMS());
        System.out.println("\t- largest:\t#" + largestClusterID + " (" + largestClusterSize + " reads)");

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
        
        System.out.println("\t- orphans:\t" + numOrphans);
        if (!checkNumAltReads) {
            emptyOrphanRecords(orphanRecords, outdir, false);
        }
        
        System.out.println("Clustered reads written in " + timer.elapsedDHMS());
        System.out.println("Clustering completed in " + timer.totalElapsedDHMS());
        System.gc();
        
        //System.out.println("before: " + NumberFormat.getInstance().format(originalNumSeq) + "\tafter: " + NumberFormat.getInstance().format(seqID));
        return counts;
    }
    
    private void emptyClusterFastaBuffer(HashMap<Integer, ArrayDeque<CompressedFastaRecord>> clusterRecords, 
            String outdir, boolean append) throws IOException {
        if (!clusterRecords.isEmpty()) {
            for (Map.Entry<Integer, ArrayDeque<CompressedFastaRecord>> e : clusterRecords.entrySet()) {
                String filePath = outdir + File.separator + e.getKey() + FASTA_EXT;
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
        String filePath = outdir + File.separator + "orphans" + FASTA_EXT;
        FastaWriter fw = new FastaWriter(filePath, append);
        for (CompressedFastaRecord f : records) {
            fw.write(f.name, f.seqbits.toString());
        }
        fw.close();
        records.clear();
    }
    
    private boolean layoutBackbones(String outFastaPath) throws IOException {
        HashMap<String, Integer> lengths = new HashMap<>(); // read id : length
        HashMap<String, String> longestAlts = new HashMap<>(); // read id : longest read id
        ArrayDeque<StrandedOverlap> dovetailRecords = new ArrayDeque<>();
        HashMap<String, Integer> artifactCutIndexes = new HashMap<>(); // read id : cut index
               
        // look for containment and dovetails
        PafReader reader = new PafReader(overlapPafInputStream);
        
        boolean checkNumAltReads = false;//minNumAltReads > 0;
//        ArrayDeque<Interval> spans = new ArrayDeque<>();
//        String prevName = null;
//        int prevLen = -1;
        ArrayDeque<String> discardReadIDs = new ArrayDeque<>();
        
        for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
            reader.next(r);
            /*
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
            */
            lengths.put(r.qName, r.qLen);
            lengths.put(r.tName, r.tLen);
            
            if (r.qName.equals(r.tName)) {
                if (cutRevCompArtifact && hasReverseComplementArtifact(r) && 
                        hasLargeOverlap(r) && hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
                    int cutIndex = getReverseComplementArtifactCutIndex(r);
                    artifactCutIndexes.put(r.qName, cutIndex);
                }
            }
            else {                
                if (hasLargeOverlap(r) && hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {                    
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
                    else if ((!longestAlts.containsKey(r.qName) ||
                            !longestAlts.containsKey(r.tName)) &&
                            isDovetailPafRecord(r)) {
                        dovetailRecords.add(pafToStrandedOverlap(r));
                    }
                }
            }
        }
        reader.close();
        /*
        if (checkNumAltReads && !spans.isEmpty() && prevName != null) {
            if (minNumAltReads > getMinCoverage(spans, prevLen, maxEdgeClip, minOverlapMatches)) {
                discardReadIDs.add(prevName);
            }
        }
        */
        System.out.println("Overlapped sequences: " + NumberFormat.getInstance().format(lengths.size()));
        if (!discardReadIDs.isEmpty()) {
            System.out.println("         - discarded: " + NumberFormat.getInstance().format(discardReadIDs.size()));
        }
        
        if (cutRevCompArtifact && !artifactCutIndexes.isEmpty()) {
            System.out.println("         - artifacts: " + NumberFormat.getInstance().format(artifactCutIndexes.size()));
        }
                
        // look for longest reads
        HashSet<String> longestSet = new HashSet<>(longestAlts.values());
        longestSet.removeAll(longestAlts.keySet());
                
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
        for (StrandedOverlap r : dovetailRecords) {                
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
                    if (r.qEnd >= lengths.get(r.qName) - maxEdgeClip && r.tEnd >= lengths.get(r.tName) - maxEdgeClip) {
                        graph.addEdge(r.tName+"+", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                        graph.addEdge(r.qName+"+", r.tName+"-", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                    }
                    else if (r.tStart <= maxEdgeClip && r.qStart <= maxEdgeClip) {
                        graph.addEdge(r.qName+"-", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                        graph.addEdge(r.tName+"-", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                    }
                }
                else {
                    if (r.qEnd >= lengths.get(r.qName) - maxEdgeClip && r.tStart <= maxEdgeClip) {
                        graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                        graph.addEdge(r.tName+"-", r.qName+"-", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                    }
                    else if (r.tEnd >= lengths.get(r.tName) - maxEdgeClip && r.qStart <= maxEdgeClip) {
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
        
        return originalNumSeq > seqID;
    }

    private class OverlapCoords {
        int qStart, qEnd, tStart, tEnd;
    }
    
    private class TargetOverlap extends OverlapCoords {
        String tName;
    }
    
    private class QueryOverlap extends OverlapCoords {
        String qName;
    }
    
    private class Overlap extends OverlapCoords {
        String qName, tName;
    }
    
    private class StrandedOverlap extends Overlap {
        boolean reverseComplemented;
    }
    
    private Overlap pafToOverlap(PafRecord r) {
        Overlap o = new Overlap();
        o.qName = r.qName;
        o.tName = r.tName;
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
        HashMap<String, Integer> lengths = new HashMap<>(); // read id : length
        HashMap<String, String> longestAlts = new HashMap<>(); // read id : longest read id
        ArrayDeque<Overlap> dovetailRecords = new ArrayDeque<>();
        HashMap<String, Integer> artifactCutIndexes = new HashMap<>(); // read id : cut index
        
        // look for containment and overlaps
        PafReader reader = new PafReader(overlapPafInputStream);
        
        boolean checkNumAltReads = false;//minNumAltReads > 0;
//        ArrayDeque<Interval> spans = new ArrayDeque<>();
//        String prevName = null;
//        int prevLen = -1;
        ArrayDeque<String> discardReadIDs = new ArrayDeque<>();
        
        for (ExtendedPafRecord r = new ExtendedPafRecord(); reader.hasNext();) {
            reader.next(r);
            /*
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
            */
            lengths.put(r.qName, r.qLen);
            lengths.put(r.tName, r.tLen);
            
            if (r.qName.equals(r.tName)) {
                if (cutRevCompArtifact && hasReverseComplementArtifact(r) && 
                        hasLargeOverlap(r) && hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
                    int cutIndex = getReverseComplementArtifactCutIndex(r);
                    artifactCutIndexes.put(r.qName, cutIndex);
                }
            }
            else {
                if (hasLargeOverlap(r) && hasGoodOverlap(r) && (!hasAlignment(r) || hasGoodAlignment(r))) {
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
                    else if ((!longestAlts.containsKey(r.qName) ||
                            !longestAlts.containsKey(r.tName)) &&
                            isStrandedDovetailPafRecord(r)) {
                        dovetailRecords.add(pafToOverlap(r));
                    }
                }
            }
        }
        reader.close();
        /*
        if (checkNumAltReads && !spans.isEmpty() && prevName != null) {
            if (minNumAltReads > getMinCoverage(spans, prevLen, maxEdgeClip, minOverlapMatches)) {
                discardReadIDs.add(prevName);
            }
        }
        */
        System.out.println("Overlapped sequences: " + NumberFormat.getInstance().format(lengths.size()));
        if (!discardReadIDs.isEmpty()) {
            System.out.println("         - discarded: " + NumberFormat.getInstance().format(discardReadIDs.size()));
        }
        
        if (cutRevCompArtifact && !artifactCutIndexes.isEmpty()) {
            System.out.println("         - artifacts: " + NumberFormat.getInstance().format(artifactCutIndexes.size()));
        }
        
        // look for longest reads
        HashSet<String> longestSet = new HashSet<>(longestAlts.values());
        longestSet.removeAll(longestAlts.keySet());
        
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
        for (Overlap r : dovetailRecords) {
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

                if (r.qEnd >= lengths.get(r.qName) - maxEdgeClip && r.tStart <= maxEdgeClip) {
                    graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                }
                else if (r.tEnd >= lengths.get(r.tName) - maxEdgeClip && r.qStart <= maxEdgeClip) {
                    graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
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
        
        return originalNumSeq > seqID;
    }
    
    public static void main(String[] args) {
        //debug
    }
}
