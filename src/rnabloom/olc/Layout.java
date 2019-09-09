/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.olc;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import org.jgrapht.Graphs;
import org.jgrapht.alg.TransitiveReduction;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.PafReader;
import rnabloom.io.PafRecord;
import static rnabloom.util.SeqUtils.bytesToString;
import static rnabloom.util.SeqUtils.reverseComplement;
import static rnabloom.util.SeqUtils.stringToBytes;


/**
 *
 * @author kmnip
 */
public class Layout {
    
    private DefaultDirectedGraph<String, OverlapEdge> graph;
    private String overlapPafPath;
    private String seqFastaPath;
    private boolean stranded;
    
    public Layout(String seqFile, String pafFile, boolean stranded) {
        this.graph = new DefaultDirectedGraph<>(OverlapEdge.class);
        this.overlapPafPath = pafFile;
        this.seqFastaPath = seqFile;
        this.stranded = stranded;
    }
    
    private static final int maxEdgeClip = 100;
    private static final float minAlnId = 0.50f;
    private static final int minOverlapMatches = 200;
    
    private class OverlapEdge extends DefaultEdge {
        public int sourceStart, sourceEnd, sinkStart, sinkEnd;
        
        public OverlapEdge(int sourceStart, int sourceEnd, int sinkStart, int sinkEnd) {
            this.sourceStart = sourceStart;
            this.sourceEnd = sourceEnd;
            this.sinkStart = sinkStart;
            this.sinkEnd = sinkEnd;
        }
    }
        
    private boolean isContainmentPafRecord(PafRecord record) {
        return record.numMatch / (float)(record.qEnd - record.qStart) >= minAlnId &&
            record.numMatch / (float)(record.tEnd - record.tStart) >= minAlnId &&
            ((record.qStart <= maxEdgeClip && record.qLen - record.qEnd <= maxEdgeClip) ||
            (record.tStart <= maxEdgeClip && record.tLen - record.tEnd <= maxEdgeClip));
    }
    
    private boolean isStrandedContainmentPafRecord(PafRecord record) {
        return !record.reverseComplemented && isContainmentPafRecord(record);
    }
    
    private boolean isDovetailPafRecord(PafRecord record) {
        if (record.numMatch >= minOverlapMatches &&
            record.numMatch / (float)(record.qEnd - record.qStart) >= minAlnId &&
            record.numMatch / (float)(record.tEnd - record.tStart) >= minAlnId) {

            if (record.reverseComplemented) {
                return (record.qEnd >= record.qLen - maxEdgeClip && record.tEnd >= record.tLen - maxEdgeClip) ||
                        (record.tStart <= maxEdgeClip && record.qStart <= maxEdgeClip);
            }
            else {
                return (record.qEnd >= record.qLen - maxEdgeClip && record.tStart <= maxEdgeClip) ||
                        (record.tEnd >= record.tLen - maxEdgeClip && record.qStart <= maxEdgeClip);
            }

        }
        
        return false;
    }
    
    private boolean isStrandedDovetailPafRecord(PafRecord record) {
        if (!record.reverseComplemented) {
            if (record.numMatch >= minOverlapMatches &&
                record.numMatch / (float)(record.qEnd - record.qStart) >= minAlnId &&
                record.numMatch / (float)(record.tEnd - record.tStart) >= minAlnId &&
                ((record.qEnd >= record.qLen - maxEdgeClip && record.tStart <= maxEdgeClip) ||
                    (record.tEnd >= record.tLen - maxEdgeClip && record.qStart <= maxEdgeClip))) {
                return true;
            }
        }
        
        return false;
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
        String name = itr.next();
        boolean reverseComplement = name.charAt(name.length()-1) == '-';
        
        byte[] bytes = sequences.get(name.substring(0, name.length()-1));
        int start = 0;
        int end = bytes.length;
        
        while (itr.hasNext()) {
            String nextName = itr.next();
            
            OverlapEdge edge = graph.getEdge(name, nextName);
            if (reverseComplement) {
                sb.append(reverseComplement(bytes, edge.sourceEnd, end));
            }
            else {
                sb.append(bytesToString(bytes, start, edge.sourceStart));
            }
            
            name = nextName;
            bytes = sequences.get(nextName.substring(0, nextName.length()-1));
            reverseComplement = nextName.charAt(name.length()-1) == '-';
            
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
    
    private void layoutBackbones(String outFastaPath) throws IOException {
        HashMap<String, Integer> lengths = new HashMap<>(); // read id -> length
        HashMap<String, String> longestAlts = new HashMap<>(); // read id -> longest read id
        ArrayDeque<PafRecord> dovetailRecords = new ArrayDeque<>();
        
        // look for containment and dovetails
        System.out.println("Parsing overlaps...");
        PafReader reader = new PafReader(overlapPafPath);
        while (reader.hasNext()) {
            PafRecord record = reader.next();
            
            String qName = record.qName;
            int qLen = record.qLen;
            
            String tName = record.tName;
            int tLen = record.tLen;
            
            lengths.put(qName, qLen);
            lengths.put(tName, tLen);
            
            if (isContainmentPafRecord(record)) {
                String shorter, longer;
                int longerLen;

                if (qLen > tLen || (qLen == tLen && qName.compareTo(tName) > 0)) {
                    shorter = tName;
                    longer = qName;
                    longerLen = qLen;
                }
                else {
                    shorter = qName;
                    longer = tName;
                    longerLen = tLen;
                }

                if (longestAlts.containsKey(shorter)) {
                    String longest = longestAlts.get(shorter);
                    if (lengths.get(longest) < longerLen) {
                        longestAlts.put(shorter, longer);
                    }
                }
                else {
                    longestAlts.put(shorter, longer);
                }
            }
            else if (!longestAlts.containsKey(qName) &&
                    !longestAlts.containsKey(tName) &&
                    isDovetailPafRecord(record)) {
                dovetailRecords.add(record);
            }
        }
        reader.close();
        
        System.out.println("Number of reads: " + lengths.size());
        System.out.println("    - dovetail:  " + dovetailRecords.size());
        
        // look for longest reads
        HashSet<String> longestSet = new HashSet<>();
        for (String r : longestAlts.values()) {
            if (longestAlts.containsKey(r)) {
                longestSet.add(longestAlts.get(r));
            }
            else {
                longestSet.add(r);
            }
        }
        
        for (String r : lengths.keySet()) {
            if (!longestSet.contains(r) && !longestAlts.containsKey(r)) {
                longestSet.add(r);
            }
        }
        
        System.out.println("    - unique:    " + longestSet.size());
        System.out.println("Constructing overlap graph...");
        
        // add nodes for reads
        for (String r : longestSet) {
            graph.addVertex(r + "+");
            graph.addVertex(r + "-");
        }        
        
        // add edges for dovetails
        for (PafRecord r : dovetailRecords) {
            if (longestSet.contains(r.qName) && longestSet.contains(r.tName)) {
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
        
        // transitive reduction
        System.out.println("Transitive reduction...");
        TransitiveReduction.INSTANCE.reduce(graph);
        
        // extract longest read sequences
        System.out.println("Extracting unique sequences...");
        HashMap<String, byte[]> longestReadSeqs = new HashMap<>(longestSet.size());
        FastaReader fr = new FastaReader(seqFastaPath);
        while (fr.hasNext()) {
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            if (longestSet.contains(name)) {
                longestReadSeqs.put(name, stringToBytes(nameSeq[1], lengths.get(name)));
            }
        }
        fr.close();
        
        // layout unambiguous paths
        System.out.println("Laying out backbones...");
        HashSet<String> visitedReadNames = new HashSet<>();
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        int seqID = 0;
        for (String n : longestSet) {
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
                String header = Integer.toString(++seqID);
                if (path.size() > 1) {
                    //print(path, ',')
                    header += " path=[" + String.join(",", path) + "]";
                }
                fw.write(header, backbone);
                
                for (String p : path) {
                    visitedReadNames.add(p.substring(0, p.length()-1));
                }
            }
        }
        fw.close();
    }
    
    private void layoutStrandedBackbones(String outFastaPath) throws IOException {
        HashMap<String, Integer> lengths = new HashMap<>(); // read id -> length
        HashMap<String, String> longestAlts = new HashMap<>(); // read id -> longest read id
        ArrayDeque<PafRecord> dovetailRecords = new ArrayDeque<>();
        
        // look for containment and overlaps
        System.out.println("Parsing overlaps...");
        PafReader reader = new PafReader(overlapPafPath);
        while (reader.hasNext()) {
            PafRecord record = reader.next();
            
            lengths.put(record.qName, record.qLen);
            lengths.put(record.tName, record.tLen);
            
            if (isStrandedContainmentPafRecord(record)) {
                String shorter, longer;
                int longerLen;
               
                if (record.qLen > record.tLen || (record.qLen == record.tLen && record.qName.compareTo(record.tName) > 0)) {
                    shorter = record.tName;
                    longer = record.qName;
                    longerLen = record.qLen;
                } 
                else {
                    shorter = record.qName;
                    longer = record.tName;
                    longerLen = record.tLen;
                }

                if (longestAlts.containsKey(shorter)) {
                    String longest = longestAlts.get(shorter);
                    if (lengths.get(longest) < longerLen) {
                        longestAlts.put(shorter, longer);
                    }
                }
                else {
                    longestAlts.put(shorter, longer);
                }
            }
            else if (!longestAlts.containsKey(record.qName) &&
                    !longestAlts.containsKey(record.tName) &&
                    isStrandedDovetailPafRecord(record)) {
                dovetailRecords.add(record);
            }
        }
        reader.close();
        
        System.out.println("Number of reads: " + lengths.size());
        System.out.println("    - dovetail:  " + dovetailRecords.size());
        
        // look for longest reads
        HashSet<String> longestSet = new HashSet<>();
        for (String r : longestAlts.values()) {
            if (longestAlts.containsKey(r)) {
                longestSet.add(longestAlts.get(r));
            }
            else {
                longestSet.add(r);
            }
        }
        
        for (String r : lengths.keySet()) {
            if (!longestSet.contains(r) && !longestAlts.containsKey(r)) {
                longestSet.add(r);
            }
        }
        
        System.out.println("    - unique:    " + longestSet.size());
        System.out.println("Constructing overlap graph...");
        
        // add nodes for reads
        for (String r : longestSet) {
            graph.addVertex(r);
        }        
        
        // add edges for dovetails
        for (PafRecord r : dovetailRecords) {
            if (!r.reverseComplemented) {
                if (longestSet.contains(r.qName) && longestSet.contains(r.tName)) {
                    if (r.qEnd >= r.qLen - maxEdgeClip && r.tStart <= maxEdgeClip) {
                        graph.addEdge(r.qName+"+", r.tName+"+", new OverlapEdge(r.qStart, r.qEnd, r.tStart, r.tEnd));
                    }
                    else if (r.tEnd >= r.tLen - maxEdgeClip && r.qStart <= maxEdgeClip) {
                        graph.addEdge(r.tName+"+", r.qName+"+", new OverlapEdge(r.tStart, r.tEnd, r.qStart, r.qEnd));
                    }
                }
            }
        }
        
        // transitive reduction
        System.out.println("Transitive reduction...");
        TransitiveReduction.INSTANCE.reduce(graph);
        
        // extract longest read sequences
        System.out.println("Extracting unique sequences...");
        HashMap<String, byte[]> longestReadSeqs = new HashMap<>(longestSet.size());
        FastaReader fr = new FastaReader(seqFastaPath);
        while (fr.hasNext()) {
            String[] nameSeq = fr.nextWithName();
            String name = nameSeq[0];
            if (longestSet.contains(name)) {
                String seq = nameSeq[1];
                longestReadSeqs.put(name, stringToBytes(seq, seq.length()));
            }
        }
        fr.close();
        
        // layout unambiguous paths
        System.out.println("Laying out backbones...");
        HashSet<String> visitedReadNames = new HashSet<>();
        FastaWriter fw = new FastaWriter(outFastaPath, false);
        int seqID = 0;
        for (String n : longestSet) {
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
                String header = Integer.toString(++seqID);
                if (path.size() > 1) {
                    //print(path, ',')
                    header += " path=[" + String.join(",", path) + "]";
                }
                fw.write(header, backbone);
                
                for (String p : path) {
                    visitedReadNames.add(p.substring(0, p.length()-1));
                }
            }
        }
        fw.close();
    }
    
    public static void main(String[] args) {
        boolean stranded = false;
        String dir = "";
        String id = "";
        String seqFastaPath = dir + "/" + id + ".fa";
        String overlapPafPath = dir + "/" + id + "_ava.paf.gz";
        String backboneFastaPath = dir + "/" + id + "_backbones.fa";
        
        try {
            Layout myLayout = new Layout(seqFastaPath, overlapPafPath, stranded);
            myLayout.writeBackboneSequences(backboneFastaPath);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
