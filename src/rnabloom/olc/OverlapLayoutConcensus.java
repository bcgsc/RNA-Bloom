/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.olc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import rnabloom.io.FastaReader;
import rnabloom.io.FastaWriter;
import rnabloom.io.PafReader;

/**
 *
 * @author kmnip
 */
public class OverlapLayoutConcensus {
    
    private final static String GZIP_EXTENSION = ".gz";
    private final static String LOG_EXTENSION = ".log";
    
    private final static String MINIMAP2 = "minimap2";
    private final static String MINIASM = "miniasm";
    private final static String RACON = "racon";
    
    private static boolean runCommand(List<String> command, String logPath) {
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);
            if (logPath != null) {
                // write stdout and stderr to file to avoid hanging due to buffer overflow
                pb.redirectErrorStream(true);
                File logFile = new File(logPath);
                pb.redirectOutput(Redirect.appendTo(logFile));
            }
            Process process = pb.start();
            int exitStatus = process.waitFor();
            return exitStatus == 0;
        }
        catch (IOException | InterruptedException e) {
            return false;
        }
    }
    
    private static boolean runCommand(List<String> command, String logPath, String gzipOutPath) {
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);
            
            if (logPath != null) {
                File logFile = new File(logPath);
                pb.redirectError(Redirect.appendTo(logFile));
            }
            
            Process process = pb.start();

            File outFile = new File(gzipOutPath);
            GZIPOutputStream zip = new GZIPOutputStream(new FileOutputStream(outFile));
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(zip, "UTF-8"));
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            
            String line;
            while ((line = reader.readLine()) != null) {
                writer.append(line);
            }
            
            int exitStatus = process.waitFor();
            
            reader.close();
            writer.close();
            
            return exitStatus == 0;
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
        
        return runCommand(command, null);
    }
    
    public static boolean hasMiniasm() {
        ArrayList<String> command = new ArrayList<>();
        command.add(MINIASM);
        command.add("-V");
        
        return runCommand(command, null);
    }
    
    public static boolean hasRacon() {
        ArrayList<String> command = new ArrayList<>();
        command.add(RACON);
        command.add("--version");
        
        return runCommand(command, null);
    }
    
    public static boolean overlapWithMinimap(String seqFastaPath, String outFastaPath, int numThreads, String options) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        command.add(MINIMAP2 + " -x ava-ont " + options + " -t " + numThreads + " " + seqFastaPath + " " + seqFastaPath + " | gzip -c > " + outFastaPath);
        
        return runCommand(command, outFastaPath + LOG_EXTENSION);
    }
    
    public static boolean mapWithMinimap(String queryFastaPath, String targetFastaPath, String outPafPath, int numThreads, String options) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");

        command.add(MINIMAP2 + " -x map-ont -c " + options + " -t " + numThreads + " " + targetFastaPath + " " + queryFastaPath + " | gzip -c > " + outPafPath);
        
        return runCommand(command, outPafPath + LOG_EXTENSION);
    }
    
    public static boolean layout(String seqFastaPath, String overlapPafPath, String backboneFastaPath, boolean stranded) {
        try {
            Layout myLayout = new Layout(seqFastaPath, overlapPafPath, stranded);
            myLayout.writeBackboneSequences(backboneFastaPath);
        } catch (Exception ex) {
            ex.printStackTrace();
            return false;
        }
        
        return true;
    }
    
    public static boolean layoutWithMiniasm(String seqFastaPath, String overlapPafPath, String layoutGfaPath) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        command.add(MINIASM + " -c 1 -e 1 -s 200 -h 100 -d 10 -g 10 -f " + seqFastaPath + " " + overlapPafPath + " | gzip -c > " + layoutGfaPath);
        
        return runCommand(command, layoutGfaPath + LOG_EXTENSION);
    }
    
    public static boolean concensusWithRacon(String queryFastaPath, String targetFastaPath, String mappingPafPath, String concensusFastaPath, int numThreads) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        command.add(RACON + " --no-trimming -u -t " + numThreads + " " + queryFastaPath + " " + mappingPafPath + " " + targetFastaPath + " > " + concensusFastaPath);
        
        return runCommand(command, concensusFastaPath + LOG_EXTENSION);
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
    
    public static boolean overlapLayoutConcensus(String readsPath, String tmpPrefix, String concensusPath, int numThreads, boolean stranded, String minimapOptions) throws IOException {
        String avaPaf = tmpPrefix + "_ava.paf.gz";
//        String backbonesGfa = tmpPrefix + "_backbones.gfa.gz";
        String backbonesFa = tmpPrefix + "_backbones.fa";
        String mapPaf = tmpPrefix + "_map.paf.gz";
        
        if (hasOnlyOneSequence(readsPath)) {
            Files.copy(Paths.get(readsPath), Paths.get(concensusPath));
            return true;
        }
        
        if (!overlapWithMinimap(readsPath, avaPaf, numThreads, minimapOptions)) {
            return false;
        }
        
        PafReader reader = new PafReader(avaPaf);
        boolean nonEmptyPafFile = reader.hasNext();
        reader.close();
        
        if (nonEmptyPafFile) {
            // lay out backbones
            if (!layout(readsPath, avaPaf, backbonesFa, stranded)) {
                return false;
            }
        }
        else {
            // PAF file is empty
            Files.copy(Paths.get(readsPath), Paths.get(concensusPath));
            return true;
        }
        
        /*
        if (!layout(readsPath, avaPaf, backbonesGfa)) {
            return false;
        }
        
        int numBackbones = convertGfaToFasta(backbonesGfa, backbonesFa);
        if (numBackbones == 0) {
            String[] nameSeq = findLongestSequence(readsPath);
            FastaWriter writer = new FastaWriter(backbonesFa, false);
            writer.write("longestread", nameSeq[1]);
            writer.close();
        }
        */
        
        if (!mapWithMinimap(readsPath, backbonesFa, mapPaf, numThreads, minimapOptions)) {
            return false;
        }
        
        return concensusWithRacon(readsPath, backbonesFa, mapPaf, concensusPath, numThreads);
    }
}
