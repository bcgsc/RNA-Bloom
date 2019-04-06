/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.olc;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;

/**
 *
 * @author kmnip
 */
public class OverlapLayoutConcensus {
    
    private final static String GZIP_EXTENSION = ".gz";
    
    private final static String MINIMAP2 = "minimap2";
    private final static String MINIASM = "miniasm";
    private final static String RACON = "racon";
    
    private static boolean runCommand(List<String> command, boolean printErrorStream) {
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);            
            Process process = pb.start();
            
            //printErrorStream = true;
            if (printErrorStream) {
                InputStream is = process.getErrorStream();
                InputStreamReader isr = new InputStreamReader(is);
                BufferedReader br = new BufferedReader(isr);
                String line;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);
                }
            }
            
            int exitVal = process.waitFor();
            return exitVal == 0;
        }
        catch (IOException | InterruptedException e) {
            return false;
        }
    }
    
    public static boolean hasOnlyOneSequence(String fasta) throws IOException {
        FastaReader reader = new FastaReader(fasta);
        int numSeq = 0;
        while (reader.hasNext()) {
            reader.next();
            if (++numSeq > 1) {
                return false;
            }
        }
        reader.close();
        
        return numSeq == 1;
    }
    
    public static boolean hasMinimap2() {
        ArrayList<String> command = new ArrayList<>();
        command.add(MINIMAP2);
        command.add("--version");
        
        return runCommand(command, false);
    }
    
    public static boolean hasMiniasm() {
        ArrayList<String> command = new ArrayList<>();
        command.add(MINIASM);
        command.add("-V");
        
        return runCommand(command, false);
    }
    
    public static boolean hasRacon() {
        ArrayList<String> command = new ArrayList<>();
        command.add(RACON);
        command.add("--version");
        
        return runCommand(command, false);
    }
    
    public static boolean overlapONT(String sequencePath, String outputPath, int numThreads) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        // CIGAR string is required for more accurate assembly with miniasm
        command.add(MINIMAP2 + " -x ava-ont -c -t " + numThreads + " " + sequencePath + " " + sequencePath + " | gzip -c > " + outputPath);
        
        return runCommand(command, false);
    }
    
    public static boolean mapONT(String queryPath, String targetPath, String outputPath, int numThreads) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        // long CIGAR string will cause a segfault in RACON
        command.add(MINIMAP2 + " -x map-ont -t " + numThreads + " " + targetPath + " " + queryPath + " | gzip -c > " + outputPath);
        
        return runCommand(command, false);
    }
    
    public static boolean layout(String sequencePath, String overlapPath, String layoutPath) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        command.add(MINIASM + " -c 1 -e 1 -s 200 -h 100 -d 10 -g 10 -f " + sequencePath + " " + overlapPath + " | gzip -c > " + layoutPath);
        
        return runCommand(command, false);
    }
    
    public static boolean concensus(String queryPath, String targetPath, String mappingPath, String concensusPath, int numThreads) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        command.add(RACON + " " + queryPath + " " + mappingPath + " " + targetPath + " -u -t " + numThreads + " > " + concensusPath);
        
        return runCommand(command, false);
    }
    
    public static int convertGfaToFasta(String gfaPath, String fastaPath) throws IOException {
        FastaWriter writer = new FastaWriter(fastaPath, false);
        
        BufferedReader br;
        if (gfaPath.endsWith(GZIP_EXTENSION) || Files.probeContentType(Paths.get(gfaPath)).equals("application/gzip")) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(gfaPath))));
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(gfaPath)));
        }
        
        int numSegments = 0;
        Iterator<String> itr = br.lines().iterator();
        while (itr.hasNext()) {
            String line = itr.next();
            if (!line.isEmpty()) {
                String[] items = line.split("\t");
                if (items[0].equals("S")) {
                    ++numSegments;
                    writer.write(items[1], items[2]);
                }
            }
        }
        
        br.close();
        writer.close();
        
        return numSegments;
    }
    
    public static String[] findLongestSequence(String fasta) throws IOException {
        FastaReader reader = new FastaReader(fasta);
        
        String[] longestNameSeq = null;
        int longestLength = -1;
        
        while (reader.hasNext()) {
            String[] nameSeq = reader.nextWithName();
            int length = nameSeq[1].length();
            if (longestLength < length) {
                longestNameSeq = nameSeq;
                longestLength = length;
            }
        }
        
        reader.close();
        
        return longestNameSeq;
    }
    
    public static boolean overlapLayoutConcensus(String readsPath, String tmpPrefix, String concensusPath, int numThreads) throws IOException {
        String avaPaf = tmpPrefix + "_ava.paf.gz";
        String gfa = tmpPrefix + "_backbones.gfa.gz";
        String gfafa = tmpPrefix + "_backbones.fa";
        String mapPaf = tmpPrefix + "_map.paf.gz";
        
        if (hasOnlyOneSequence(readsPath)) {
            Files.copy(Paths.get(readsPath), Paths.get(gfafa));
            return true;
        }
        
        if (!overlapONT(readsPath, avaPaf, numThreads)) {
            return false;
        }
        
        if (!layout(readsPath, avaPaf, gfa)) {
            return false;
        }
        
        int numBackbones = convertGfaToFasta(gfa, gfafa);
        if (numBackbones == 0) {
            String[] nameSeq = findLongestSequence(readsPath);
            FastaWriter writer = new FastaWriter(gfafa, false);
            writer.write(nameSeq[0], nameSeq[1]);
            writer.close();
        }
        
        if (!mapONT(readsPath, gfafa, mapPaf, numThreads)) {
            return false;
        }
        
        return concensus(readsPath, gfafa, mapPaf, concensusPath, numThreads);
    }
}
