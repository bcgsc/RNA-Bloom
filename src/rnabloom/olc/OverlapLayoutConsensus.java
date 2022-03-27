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
import java.io.Writer;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Files;
import java.util.ArrayDeque;
import java.util.ArrayList;
import static rnabloom.io.Constants.FASTA_EXT;
import rnabloom.io.ExtendedPafRecord;
import rnabloom.io.PafReader;
import static rnabloom.util.CommandLine.runCommand;
import static rnabloom.util.Common.convertToRoundedPercent;
import static rnabloom.util.FileUtils.deleteIfExists;
import static rnabloom.util.FileUtils.hasOnlyOneSequence;
import static rnabloom.util.FileUtils.readIntArrayFromFile;
import static rnabloom.util.FileUtils.symlinkRemoveExisting;
import static rnabloom.util.FileUtils.touch;
import static rnabloom.util.FileUtils.writeIntArrayToFile;
import static rnabloom.util.PafUtils.hasGoodAlignment;
import static rnabloom.util.PafUtils.hasLargeOverlap;
import static rnabloom.util.PafUtils.isContainmentPafRecord;
import rnabloom.util.PolyATailFinder;
import rnabloom.util.Timer;
import static rnabloom.util.FileUtils.getTextFileWriter;

/**
 *
 * @author Ka Ming Nip
 */
public class OverlapLayoutConsensus {
    
    private final static String GZIP_EXTENSION = ".gz";
    private final static String LOG_EXTENSION = ".log";
    
    private final static String MINIMAP2 = "minimap2";
    private final static String RACON = "racon";
    
    private final static String PRESET_PACBIO = "pb";
    private final static String PRESET_ONT = "ont";
    
    public static enum STATUS {SUCCESS, EMPTY, FAIL};
                
    public static boolean hasMinimap2() {
        ArrayList<String> command = new ArrayList<>();
        command.add(MINIMAP2);
        command.add("--version");
        
        return runCommand(command, null);
    }
        
    public static boolean hasRacon() {
        ArrayList<String> command = new ArrayList<>();
        command.add(RACON);
        command.add("--version");
        
        return runCommand(command, null);
    }
    
    public static boolean overlapWithMinimap(String seqFastaPath, String outPafPath, int numThreads, 
            boolean align, String minimapOptions, boolean usePacBioPreset) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
                
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
                
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (!minimapOptions.contains("-X")) {
            minimapOptions += " -X"; 
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x ava-" + preset + " " + minimapOptions + " " + seqFastaPath + " " + seqFastaPath + " | gzip --fast -c > " + outPafPath);
        
        return runCommand(command, outPafPath + LOG_EXTENSION);
    }
    
    public static STATUS overlapWithMinimapAndExtractUnique(String seqFastaPath, String uniqueFastaPath,
            int numThreads, boolean align, String minimapOptions, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            int minSeqDepth, boolean usePacBioPreset, boolean verbose, PolyATailFinder polyaFinder) {
        if (verbose) {
            System.out.println("Overlapping sequences...");
        }
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
                
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (!minimapOptions.contains("-X")) {
            minimapOptions += " -X"; 
        }
        
        if (stranded && !minimapOptions.contains("--for-only")) {
            minimapOptions += " --for-only";
        }
                
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x ava-" + preset + " " + minimapOptions + " " + seqFastaPath + " " + seqFastaPath);
        
        STATUS status = STATUS.FAIL;
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(uniqueFastaPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            Layout myLayout = new Layout(seqFastaPath, process.getInputStream(), stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, false, minSeqDepth, verbose, polyaFinder);
            long numUnique = myLayout.extractUniqueFromOverlaps(uniqueFastaPath);
            if (numUnique <= 0) {
                return STATUS.EMPTY;
            }
            
            int exitStatus = process.waitFor();
            if (exitStatus == 0) {
                status = STATUS.SUCCESS;
            }
            else {
                status = STATUS.FAIL;
            }
        }
        catch (IOException | InterruptedException e) {
            status = STATUS.FAIL;
        }
        
        return status;
    }
    
    public static boolean trimSplitByReadDepth(String queryFastaPath, 
            String targetFastaPath, String outFastaPath,
            int numThreads, boolean align, String minimapOptions, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth, boolean usePacBioPreset, boolean verbose, int minLen) {
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
                
//        if (!minimapOptions.contains("-g")) {
//            minimapOptions += " -g 300";
//        }
//        
//        if (!minimapOptions.contains("-r")) {
//            minimapOptions += " -r " + maxIndelSize;
//        }
        
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (stranded && !minimapOptions.contains("--for-only")) {
            minimapOptions += " --for-only";
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x map-" + preset + " " + minimapOptions + " " + targetFastaPath + " " + queryFastaPath);
        
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(outFastaPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();
            
            Layout myLayout = new Layout(queryFastaPath, process.getInputStream(), stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose, null);
            myLayout.trimSplitByReadDepth(targetFastaPath, outFastaPath, minLen);
            
            int exitStatus = process.waitFor();
            if (exitStatus != 0) {
                return false;
            }    
        }
        catch (IOException | InterruptedException e) {
            return false;
        }
        
        return true;
    }
    
    public static boolean mapWithMinimapAndExtractUnique(String queryFastaPath,
            String targetFastaPath, String uniqueFastaPath,
            int numThreads, boolean align, String minimapOptions, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth, boolean usePacBioPreset,
            boolean verbose, PolyATailFinder polyaFinder) {
        if (verbose) {
            System.out.println("Mapping sequences...");
        }
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
                
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (stranded && !minimapOptions.contains("--for-only")) {
            minimapOptions += " --for-only";
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x map-" + preset + " " + minimapOptions + " " + targetFastaPath + " " + queryFastaPath);
        
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(uniqueFastaPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();
            
            Layout myLayout = new Layout(queryFastaPath, process.getInputStream(), stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose, polyaFinder);
            myLayout.extractUniqueFromMapping(uniqueFastaPath);
            
            int exitStatus = process.waitFor();
            if (exitStatus != 0) {
               return false;
            }
        }
        catch (IOException | InterruptedException e) {
            return false;
        }
        
        return true;
    }
    
    public static int[] mapWithMinimapAndExtractClusters(String queryFastaPath, 
            String targetFastaPath, String clustersDir,
            int numThreads, boolean align, String minimapOptions, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth, boolean usePacBioPreset, boolean verbose) {
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
        
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (stranded && !minimapOptions.contains("--for-only")) {
            minimapOptions += " --for-only";
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x map-" + preset + " " + minimapOptions + " " + targetFastaPath + " " + queryFastaPath);
        
        int[] clusterSizes = null;
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(clustersDir + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();
            
            Layout myLayout = new Layout(queryFastaPath, process.getInputStream(), stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose, null);
            clusterSizes = myLayout.extractClustersFromMapping(clustersDir);
            
            int exitStatus = process.waitFor();
            if (exitStatus != 0) {
                clusterSizes = null;
            }    
        }
        catch (IOException | InterruptedException e) {
            return null;
        }
        
        return clusterSizes;
    }
    
    public static int[] overlapWithMinimapAndExtractClusters(String seqFastaPath, String clusterdir,
            int numThreads, boolean align, String minimapOptions, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, int maxMergedClusterSize, boolean verbose) {
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
                
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }

//        if (!minimapOptions.contains("-g")) {
//            minimapOptions += " -g 300";
//        }
        
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (!minimapOptions.contains("-X")) {
            minimapOptions += " -X"; 
        }
        
        if (stranded && !minimapOptions.contains("--for-only")) {
            minimapOptions += " --for-only";
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x ava-" + preset + " " + minimapOptions + " " + seqFastaPath + " " + seqFastaPath);
        
        int[] clusterSizes = null;
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(clusterdir + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();
            
            Layout myLayout = new Layout(seqFastaPath, process.getInputStream(), stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose, null);
            clusterSizes = myLayout.extractClustersFromOverlaps(clusterdir, maxMergedClusterSize);
            
            int exitStatus = process.waitFor();
            if (exitStatus != 0) {
                clusterSizes = null;
            }    
        }
        catch (IOException | InterruptedException e) {
            return null;
        }
        
        return clusterSizes;
    }
    
    public static STATUS overlapWithMinimapAndLayout(String seqFastaPath, String layoutFastaPath,
            int numThreads, boolean align, String minimapOptions, boolean stranded, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, boolean verbose) {
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
                
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
        
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (!cutRevCompArtifact) {
            if (!minimapOptions.contains("-X")) {
                minimapOptions += " -X"; 
            }

            if (stranded && !minimapOptions.contains("--for-only")) {
                minimapOptions += " --for-only";
            }
        }
        
//        if (maxIndelSize >= 0 && !minimapOptions.contains("-r ")) {
//            minimapOptions += " -r " + maxIndelSize;
//        }
        
//        if (!minimapOptions.contains("-f ")) {
//            minimapOptions += " -f 0.0001";
//        }
                
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x ava-" + preset + " " + minimapOptions + " " + seqFastaPath + " " + seqFastaPath);
        
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(layoutFastaPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            boolean reduced = layout(seqFastaPath, process.getInputStream(), layoutFastaPath, stranded, maxEdgeClip, 
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose);
            
            int exitStatus = process.waitFor();
            
            if (exitStatus != 0) {
                return STATUS.FAIL;
            }
            else if (!reduced) {
                return STATUS.EMPTY;
            }
            else {
                return STATUS.SUCCESS;
            }
        }
        catch (IOException | InterruptedException e) {
            return STATUS.FAIL;
        }
    }

    public static STATUS overlapWithMinimapAndLayoutSimple(String seqFastaPath, String layoutFastaPath,
            int numThreads, boolean align, String minimapOptions, boolean stranded, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, boolean verbose, PolyATailFinder polyaFinder) {
        if (verbose) {
            System.out.println("Overlapping sequences...");
        }
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
                
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
        
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (!minimapOptions.contains("-X")) {
            minimapOptions += " -X"; 
        }
        
        if (stranded && !minimapOptions.contains("--for-only")) {
            minimapOptions += " --for-only";
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x ava-" + preset + " " + minimapOptions + " " + seqFastaPath + " " + seqFastaPath);
        
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(layoutFastaPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            boolean reduced = layoutSimple(seqFastaPath, process.getInputStream(), layoutFastaPath, stranded, maxEdgeClip, 
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose,
                    polyaFinder);
            
            int exitStatus = process.waitFor();
            
            if (exitStatus != 0) {
                return STATUS.FAIL;
            }
            else if (!reduced) {
                return STATUS.EMPTY;
            }
            else {
                return STATUS.SUCCESS;
            }
        }
        catch (IOException | InterruptedException e) {
            return STATUS.FAIL;
        }
    }
    
    public static STATUS overlapWithMinimapAndLayoutGreedy(String seqFastaPath, String layoutFastaPath,
            int numThreads, boolean align, String minimapOptions, boolean stranded, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, String mappingPafPath, String polyAReadNamesPath,
            boolean verbose, PolyATailFinder polyaFinder, String sampleReadLengthsPath,
            String txptNamePrefix, boolean writeUracil) {
        if (verbose) {
            System.out.println("Overlapping sequences...");
        }
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align && !minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
                
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (!minimapOptions.contains("-X")) {
            minimapOptions += " -X"; 
        }
        
        if (stranded && !minimapOptions.contains("--for-only")) {
            minimapOptions += " --for-only";
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x ava-" + preset + " " + minimapOptions + " " + seqFastaPath + " " + seqFastaPath);
        
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(layoutFastaPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            boolean reduced = layoutGreedy(seqFastaPath, process.getInputStream(), layoutFastaPath, stranded, maxEdgeClip, 
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, mappingPafPath,
                    polyAReadNamesPath, verbose, polyaFinder, sampleReadLengthsPath, txptNamePrefix,
                    writeUracil);
            
            int exitStatus = process.waitFor();
            
            if (exitStatus != 0) {
                return STATUS.FAIL;
            }
            else if (!reduced) {
                return STATUS.EMPTY;
            }
            else {
                return STATUS.SUCCESS;
            }
        }
        catch (IOException | InterruptedException e) {
            return STATUS.FAIL;
        }
    }
    
    public static boolean mapWithMinimap(String queryFastaPath, String targetFastaPath, String outPafPath,
            int numThreads, String minimapOptions, boolean usePacBioPreset) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");

        if (!minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
        
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x map-" + preset + " " + minimapOptions + " " + targetFastaPath + " " + queryFastaPath + " | gzip --fast -c > " + outPafPath);
        
        return runCommand(command, outPafPath + LOG_EXTENSION);
    }
    
    public static STATUS mapWithMinimapFiltered(String queryFastaPath, String targetFastaPath, String outPafPath,
            int numThreads, String minimapOptions, boolean usePacBioPreset, boolean stranded, int maxIndelSize, 
            int minOverlapMatches, float minAlnId, int maxEdgeClip, boolean verbose) {
        Timer timer = new Timer();
        if (verbose) {
            System.out.println("Mapping sequences...");
        }
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");

        if (!minimapOptions.contains("-c")) {
            minimapOptions += " -c";
        }
        
        if (!minimapOptions.contains("-N")) {
            minimapOptions += " -N 10";
        }
        
        if (!minimapOptions.contains("-p")) {
            minimapOptions += " -p 0.1";
        }
        
        if (numThreads > 0) {
            if (numThreads >= 2) {
                --numThreads;
            }
            
            if (!minimapOptions.contains("-t")) {
                minimapOptions += " -t " + numThreads;
            }
        }
        
        if (stranded && !minimapOptions.contains("--for-only")) {
            minimapOptions += " --for-only";
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x map-" + preset + " " + minimapOptions + " " + targetFastaPath + " " + queryFastaPath);

        try {
            Writer bw = getTextFileWriter(outPafPath, false);
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(outPafPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            long numGoodRecords = 0;
            PafReader reader = new PafReader(process.getInputStream());
            String qName = null;
            ArrayDeque<ExtendedPafRecord> currentRecords = new ArrayDeque<>();
            ExtendedPafRecord primary = null;
            for (ExtendedPafRecord record; reader.hasNext(); ) {
                record = reader.next();
                if (!record.qName.equals(qName)) {
                    qName = record.qName;
                    
                    if (!currentRecords.isEmpty()) {
                        if (primary != null) {
                            float id = primary.numMatch/(float)primary.blockLen;
                            // filter secondary alignments based on primary alignment's identity
                            float threshold = id * id;
                            for (ExtendedPafRecord r : currentRecords) {
                                if (r.numMatch/(float)r.blockLen >= threshold) {
                                    bw.write(r.toString());
                                    bw.write('\n');
                                }
                            }
                        }
                        else {
                            for (ExtendedPafRecord r : currentRecords) {
                                bw.write(r.toString());
                                bw.write('\n');
                            }
                        }
                    }
                    
                    currentRecords = new ArrayDeque<>();
                    primary = null;
                }
                
                if (!stranded || !record.reverseComplemented) {
                    if ((hasLargeOverlap(record, minOverlapMatches) ||
                            isContainmentPafRecord(record, maxEdgeClip)) &&
                            hasGoodAlignment(record, maxIndelSize, minAlnId)) {
                        currentRecords.add(record);
                        ++numGoodRecords;
                    }
                }
                
                if (record.isPrimary) {
                    primary = record;
                }
            }
            
            // process the final batch
            if (!currentRecords.isEmpty()) {
                if (primary != null) {
                    float id = primary.numMatch/(float)primary.blockLen;
                    float threshold = id * id;
                    for (ExtendedPafRecord r : currentRecords) {
                        if (r.numMatch/(float)r.blockLen >= threshold) {
                            bw.write(r.toString());
                            bw.write('\n');
                        }
                    }
                }
                else {
                    for (ExtendedPafRecord r : currentRecords) {
                        bw.write(r.toString());
                        bw.write('\n');
                    }
                }
            }
            
            reader.close();
            bw.close();
            
            int exitStatus = process.waitFor();
            if (exitStatus != 0) {
                return STATUS.FAIL;
            }
            else if (numGoodRecords == 0) {
                return STATUS.EMPTY;
            }
        }
        catch (IOException | InterruptedException e) {
            return STATUS.FAIL;
        }
        
        if (verbose) {
            System.out.println("Mapping completed in " + timer.elapsedDHMS());
        }
        
        return STATUS.SUCCESS;
    }
        
    public static boolean layout(String seqFastaPath, InputStream overlapPafInputStream, String backboneFastaPath,
            boolean stranded, int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth, boolean verbose) {
        try {
            Layout myLayout = new Layout(seqFastaPath, overlapPafInputStream, stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose, null);
            return myLayout.layoutBackbones(backboneFastaPath);
        } catch (Exception ex) {
            ex.printStackTrace();
            return false;
        }
    }

    public static boolean layoutSimple(String seqFastaPath, InputStream overlapPafInputStream, String backboneFastaPath,
            boolean stranded, int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth, boolean verbose, PolyATailFinder polyaFinder) {
        try {
            Layout myLayout = new Layout(seqFastaPath, overlapPafInputStream, stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose, polyaFinder);
            myLayout.extractSimplePaths(backboneFastaPath);
            return true;
        } catch (Exception ex) {
            ex.printStackTrace();
            return false;
        }
    }
    
    public static boolean layoutGreedy(String seqFastaPath, InputStream overlapPafInputStream, String backboneFastaPath,
            boolean stranded, int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth, String mappingPafPath, String polyAReadNamesPath,
            boolean verbose, PolyATailFinder polyaFinder, String sampleReadLengthsPath,
            String txptNamePrefix, boolean writeUracil) {
        try {
            Layout myLayout = new Layout(seqFastaPath, overlapPafInputStream, stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth, verbose, polyaFinder);
            myLayout.extractGreedyPaths(backboneFastaPath, mappingPafPath,
                    polyAReadNamesPath, sampleReadLengthsPath, txptNamePrefix, writeUracil);
            return true;
        } catch (Exception ex) {
            ex.printStackTrace();
            return false;
        }
    }
        
    public static boolean consensusWithRacon(String queryFastaPath, String targetFastaPath, 
            String mappingPafPath, String consensusFastaPath, int numThreads, 
            boolean keepUnpolished, boolean verbose) {
        Timer timer = new Timer();
        if (verbose) {
            System.out.println("Polishing sequences...");
        }
        
        String extraOptions = "";
        if (keepUnpolished) {
            extraOptions += " --no-trimming -u ";
        }
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        command.add(RACON + extraOptions + " -t " + numThreads + " " + queryFastaPath + " " + mappingPafPath + " " + targetFastaPath + " > " + consensusFastaPath);
        //-f
        //--no-trimming -u
        
        boolean status = runCommand(command, consensusFastaPath + LOG_EXTENSION);
        
        if (verbose) {
            System.out.println("Polishing completed in " + timer.elapsedDHMS());
        }
        
        return status;
    }
        
    public static boolean overlapLayout(String readsPath, String layoutPath,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, boolean verbose) throws IOException {
        
        if (hasOnlyOneSequence(readsPath)) {
            symlinkRemoveExisting(readsPath, layoutPath);
            return true;
        }
        
        STATUS s = overlapWithMinimapAndLayout(readsPath, layoutPath,
            numThreads, true, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
            minSeqDepth, usePacBioPreset, verbose);
        
        switch (s) {
            case FAIL:
                return false;
            case EMPTY:
                symlinkRemoveExisting(readsPath, layoutPath);
        }
        
        return true;
    }

    public static boolean mapAndConsensus(String readsPath, String targetPath, String tmpPrefix, String consensusPath, 
            int numThreads, String minimapOptions, boolean usePacBioPreset, boolean keepUnpolished, boolean verbose) throws IOException {
        String mapPaf = tmpPrefix + "_map.paf.gz";
        deleteIfExists(mapPaf);
        
        if (!mapWithMinimap(readsPath, targetPath, mapPaf, numThreads, minimapOptions, usePacBioPreset)) {
            return false;
        }
        
        return consensusWithRacon(readsPath, targetPath, mapPaf, consensusPath, numThreads, keepUnpolished, verbose);
    }
    
    public static boolean overlapLayoutConsensus(String readsPath, String tmpPrefix, String consensusPath, 
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, boolean align, boolean keepUnpolished, boolean verbose) throws IOException {
        String backbonesFa = tmpPrefix + "_backbones.fa";
        String mapPaf = tmpPrefix + "_map.paf.gz";
        
        deleteIfExists(backbonesFa);
        deleteIfExists(mapPaf);
        
        STATUS s = overlapWithMinimapAndLayout(readsPath, backbonesFa,
            numThreads, align, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
            minSeqDepth, usePacBioPreset, verbose);
        
        if (s == STATUS.EMPTY) {
            // PAF is empty
            symlinkRemoveExisting(readsPath, consensusPath);
            return true;
        }
        
        if (!mapWithMinimap(readsPath, backbonesFa, mapPaf, numThreads,
                minimapOptions, usePacBioPreset)) {
            return false;
        }
        
//        s = mapWithMinimapFiltered(readsPath, backbonesFa, mapPaf, numThreads,
//                minimapOptions, usePacBioPreset, stranded, maxIndelSize, 
//                minOverlapMatches, minAlnId);
//        switch(s) {
//            case FAIL:
//                return false;
//            case EMPTY:
//                symlinkRemoveExisting(readsPath, consensusPath);
//                return true;
//        }
        
        return consensusWithRacon(readsPath, backbonesFa, mapPaf, consensusPath,
                numThreads, keepUnpolished, verbose);
    }
    
    public static boolean overlapLayoutConsensus2(String readsPath, String tmpPrefix, String consensusPath, 
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean usePacBioPreset, boolean alignReads, boolean alignBackbones,
            boolean keepUnpolished, boolean verbose) throws IOException {
        String backbonesFa = tmpPrefix + "_backbones.fa";
        //String polishedFa = tmpPrefix + "_polished.fa";
        String backbonesFa2 = tmpPrefix + "_backbones2.fa";
        String mapPaf = tmpPrefix + "_map.paf.gz";

        deleteIfExists(backbonesFa);
        //deleteIfExists(polishedFa);
        deleteIfExists(backbonesFa2);
        deleteIfExists(mapPaf);

//        if (hasOnlyOneSequence(readsPath)) {
//            symlinkRemoveExisting(readsPath, consensusPath);
//            return true;
//        }
        
        STATUS s = overlapWithMinimapAndLayout(readsPath, backbonesFa,
            numThreads, alignReads, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, false,
            1, usePacBioPreset, verbose);

        if (s == STATUS.EMPTY) {
            // either PAF is empty or no backbones can be made
            symlinkRemoveExisting(readsPath, consensusPath);
            return true;
        }

        if (hasOnlyOneSequence(backbonesFa)) {
            // polish backbone #1
            
            boolean status = mapWithMinimap(readsPath, backbonesFa, mapPaf, numThreads, minimapOptions, usePacBioPreset);
            if (!status) {
                return false;
            }

            return consensusWithRacon(readsPath, backbonesFa, mapPaf, consensusPath,
                    numThreads, keepUnpolished, verbose);
        }
        else {
            // layout backbone #2

            s = overlapWithMinimapAndLayout(backbonesFa, backbonesFa2,
                        numThreads, alignBackbones, minimapOptions, stranded, maxEdgeClip,
                        minAlnId, minOverlapMatches, maxIndelSize, false,
                        1, usePacBioPreset, verbose);

            if (s == STATUS.EMPTY) {
                // either PAF is empty or no backbones can be made
//                symlinkRemoveExisting(backbonesFa2, consensusPath);
//                return true;
                backbonesFa2 = backbonesFa;
            }

            // polish backbone #2

            boolean status = mapWithMinimap(readsPath, backbonesFa2, mapPaf, numThreads, minimapOptions, usePacBioPreset);
            if (!status) {
                return false;
            }

            return consensusWithRacon(readsPath, backbonesFa2, mapPaf, consensusPath,
                    numThreads, keepUnpolished, verbose);
        }
    }
    
    public static boolean seededUniqueOLC(String readsPath, String seedsPath, String outFastaPath, String tmpPrefix,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize,
            int minSeqDepth, boolean usePacBioPreset, boolean verbose,
            String polyAReadNamesPath) throws IOException {

        PolyATailFinder polyaFinder = null;
        
        String overlappedSeedsPath = tmpPrefix + "0.seeds" + FASTA_EXT + GZIP_EXTENSION;
        String uniqueFastaPath = tmpPrefix + "1.nr" + FASTA_EXT + GZIP_EXTENSION;
        String cutFastaPath = tmpPrefix + "2.cut" + FASTA_EXT + GZIP_EXTENSION;
        String unitigsFastaPath = tmpPrefix + "3.uni" + FASTA_EXT + GZIP_EXTENSION;
        String readsToSimplePafPath = tmpPrefix + "4.map.paf.gz";
        String polishedFastaPath = tmpPrefix + "5.pol" + FASTA_EXT;
        
        deleteIfExists(overlappedSeedsPath);
        deleteIfExists(uniqueFastaPath);
        deleteIfExists(cutFastaPath);
        deleteIfExists(unitigsFastaPath);
        deleteIfExists(readsToSimplePafPath);
        deleteIfExists(polishedFastaPath);
        deleteIfExists(outFastaPath);
        
//        String minimapOptionsNoGaps = minimapOptions;
//        if (!minimapOptionsNoGaps.contains("-g ")) {
//            minimapOptionsNoGaps += " -g 300";
//        }
//        if (!minimapOptionsNoGaps.contains("-r ")) {
//            minimapOptionsNoGaps += " -r " + maxIndelSize;
//        }

        // 0. overlap seeds
        boolean success = overlapLayout(seedsPath, overlappedSeedsPath,
            numThreads, stranded, minimapOptions, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, false,
            minSeqDepth, usePacBioPreset, true);
        if (!success) {
            return false;
        }
        
        System.gc();
        
        // 1. map all reads and extract unique reads
        mapWithMinimapAndExtractUnique(readsPath,
            overlappedSeedsPath, uniqueFastaPath,
            numThreads, true, minimapOptions, stranded,
            maxEdgeClip, minAlnId, minOverlapMatches, maxIndelSize,
            false, minSeqDepth, usePacBioPreset, verbose,
            polyaFinder);
        
        System.gc();
        
        // 2. overlap all reads and cut reads
        STATUS status = overlapWithMinimapAndExtractUnique(uniqueFastaPath, cutFastaPath,
            numThreads, false, minimapOptions, stranded,
            maxEdgeClip, minAlnId*minAlnId, minOverlapMatches, maxIndelSize,
            minSeqDepth, usePacBioPreset, verbose, polyaFinder);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        System.gc();
        
        // 3. overlap (with alignment) cutted reads and extract unitigs
        status = overlapWithMinimapAndLayoutSimple(cutFastaPath, unitigsFastaPath,
            numThreads, true, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, false,
            1, usePacBioPreset, verbose, polyaFinder);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        System.gc();
        
        // 4. align all reads to unitigs
        status = mapWithMinimapFiltered(readsPath, unitigsFastaPath, readsToSimplePafPath,
            numThreads, minimapOptions, usePacBioPreset, stranded, maxIndelSize,
            minOverlapMatches, minAlnId, maxEdgeClip, verbose);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        System.gc();
        
        // 5. polish unitigs
        success = consensusWithRacon(readsPath, unitigsFastaPath, 
            readsToSimplePafPath, polishedFastaPath, numThreads, true, verbose);
        if (!success) {
            return false;
        }
        
        // 6. overlap (with alignment) polished unitigs and lay out paths
        status = overlapWithMinimapAndLayoutGreedy(polishedFastaPath, outFastaPath,
            numThreads, true, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, false,
            1, usePacBioPreset, readsToSimplePafPath, polyAReadNamesPath,
            verbose, polyaFinder, null, "rnabloom", false);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        return true;
    }
    
    public static boolean uniqueOLC(String readsPath, String inFastaPath,
            String outFastaPath, String tmpPrefix, String sampleReadLengthsPath,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize,
            int minSeqDepth, boolean usePacBioPreset, boolean verbose,
            String polyAReadNamesPath, PolyATailFinder polyaFinder, String txptNamePrefix,
            boolean writeUracil) throws IOException {

        String uniqueFastaPath = tmpPrefix + "1.nr" + FASTA_EXT + GZIP_EXTENSION;
        String simpleFastaPath = tmpPrefix + "2.uni" + FASTA_EXT + GZIP_EXTENSION;
        String readsToSimplePafPath = tmpPrefix + "3.map.paf.gz";
        String polishedSimpleFastaPath = tmpPrefix + "4.pol" + FASTA_EXT;
        
        deleteIfExists(uniqueFastaPath);
        deleteIfExists(simpleFastaPath);
        deleteIfExists(readsToSimplePafPath);
        deleteIfExists(polishedSimpleFastaPath);
        deleteIfExists(outFastaPath);

        String minimapOptionsNoGaps = minimapOptions;
        if (!minimapOptionsNoGaps.contains("-g ")) {
            minimapOptionsNoGaps += " -g 300";
        }
        if (maxIndelSize >= 0 && !minimapOptionsNoGaps.contains("-r ")) {
            minimapOptionsNoGaps += " -r " + maxIndelSize;
        }
        
        System.gc();
        
        // 1. overlap all reads and extract unique reads
        STATUS status = overlapWithMinimapAndExtractUnique(inFastaPath, uniqueFastaPath,
            numThreads, false, minimapOptionsNoGaps, stranded,
            maxEdgeClip, minAlnId*minAlnId, minOverlapMatches, maxIndelSize,
            minSeqDepth, usePacBioPreset, verbose, polyaFinder);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        System.gc();
                
        // 2. overlap (with alignment) unique reads and extract unitigs
        status = overlapWithMinimapAndLayoutSimple(uniqueFastaPath, simpleFastaPath,
            numThreads, true, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, false,
            1, usePacBioPreset, verbose, polyaFinder);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        System.gc();
        
        // 3. map all reads to unitigs
        status = mapWithMinimapFiltered(readsPath, simpleFastaPath, readsToSimplePafPath,
            numThreads, minimapOptions, usePacBioPreset, stranded, maxIndelSize,
            minOverlapMatches, minAlnId, maxEdgeClip, verbose);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        System.gc();
        
        // 4. polish unitigs
        boolean keepUnpolished = true;
        boolean success = consensusWithRacon(readsPath, simpleFastaPath, 
            readsToSimplePafPath, polishedSimpleFastaPath, numThreads, keepUnpolished, verbose);
        if (!success) {
            return false;
        }
        
//        if (!minimapOptionsNoGaps.contains("-f ")) {
//            minimapOptionsNoGaps += " -f 0.0001";
//        }
        
        // 5. overlap (with alignment) unitigs and lay out paths
        status = overlapWithMinimapAndLayoutGreedy(polishedSimpleFastaPath, outFastaPath,
            numThreads, true, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, false,
            1, usePacBioPreset, readsToSimplePafPath, polyAReadNamesPath,
            verbose, polyaFinder, sampleReadLengthsPath, txptNamePrefix,
            writeUracil);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        return true;
    }
    
    public static int avaClusteredOLC(String readsPath, String clustersdir,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, int maxMergedClusterSize, boolean forceOverwrite,
            boolean verbose) throws IOException {
        
        File clusteringSizesFile = new File(clustersdir + File.separator + "cluster_sizes.txt");
        int[] clusterSizes = null;
        
        if (forceOverwrite) {
            Files.deleteIfExists(clusteringSizesFile.toPath());
            
            clusterSizes = overlapWithMinimapAndExtractClusters(readsPath, clustersdir,
                    numThreads, false, minimapOptions, stranded, maxEdgeClip,
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
                    minSeqDepth, usePacBioPreset, maxMergedClusterSize, verbose);
            
            if (clusterSizes == null || clusterSizes.length == 0) {
                return -1;
            }
            
            writeIntArrayToFile(clusteringSizesFile, clusterSizes);
        }
        else {
            if (clusteringSizesFile.exists()) {
                clusterSizes = readIntArrayFromFile(clusteringSizesFile);
                
                if (clusterSizes == null || clusterSizes.length == 0) {
                    return -1;
                }
            }
            else {
                clusterSizes = overlapWithMinimapAndExtractClusters(readsPath, clustersdir,
                    numThreads, false, minimapOptions, stranded, maxEdgeClip,
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
                    minSeqDepth, usePacBioPreset, maxMergedClusterSize, verbose);
                
                if (clusterSizes == null || clusterSizes.length == 0) {
                    return -1;
                }
                
                writeIntArrayToFile(clusteringSizesFile, clusterSizes);
            }
        }
                
        minSeqDepth = 1;
        cutRevCompArtifact = false;
//        boolean keepUnpolished = true;
        
        int numClusters = clusterSizes.length;
        for (int c=0; c<numClusters; ++c) {            
            int cid = c+1;
            int numReads = clusterSizes[c];
            
            File assemblyDoneStamp = new File(clustersdir + File.separator + cid + ".DONE");
            System.out.println("Processing " + cid + " (" + numReads + " reads)...");
            
            if (forceOverwrite || !assemblyDoneStamp.exists()) {
                String clusterPrefix = clustersdir + File.separator + cid;
                String inFastaPath = clusterPrefix + FASTA_EXT;
//                String tmpPrefix = clusterPrefix + "_tmp";
                String finalFastaPath = clusterPrefix + "_transcripts" + FASTA_EXT;

                // do not align during AVA overlap if the cluster has to many reads
                boolean align = numReads <= 500000;
                boolean status = false;
//                if (numReads >= 3) {
//                    status = overlapLayoutConsensus(inFastaPath, tmpPrefix, finalFastaPath, 
//                                    numThreads, stranded, minimapOptions, maxEdgeClip,
//                                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth,
//                                    usePacBioPreset, align, keepUnpolished);
//                }
//                else {
                    STATUS s = overlapWithMinimapAndLayout(inFastaPath, finalFastaPath, 
                                    numThreads, align, minimapOptions, stranded, maxEdgeClip,
                                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth,
                                    usePacBioPreset, verbose);
                    switch (s) {
                        case SUCCESS:
                            status = true;
                            break;
                        case FAIL:
                            status = false;
                            break;
                        case EMPTY:
                            status = true;
                            symlinkRemoveExisting(inFastaPath, finalFastaPath);
                            break;
                    }
//                }
                
                if (status) {
                    touch(assemblyDoneStamp);
                }
                else {
                    throw new IOException("Error processing cluster " + cid);
                }
            }
            else {
                System.out.println("WARNING: cluster #" + cid + " was already processed!");
            }
        }
        
        return numClusters;
    }
    
    public static int mapClusteredOLC(String readsPath, String targetPath, String clustersdir,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, boolean forceOverwrite,
            String polyAReadNamesPath, PolyATailFinder polyaFinder) throws IOException {
        
        File clusteringSizesFile = new File(clustersdir + File.separator + "cluster_sizes.txt");
        int[] clusterSizes = null;
        
        if (forceOverwrite) {
            Files.deleteIfExists(clusteringSizesFile.toPath());
            
            clusterSizes = mapWithMinimapAndExtractClusters(readsPath, targetPath, clustersdir,
                    numThreads, true, minimapOptions, stranded, maxEdgeClip,
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
                    minSeqDepth, usePacBioPreset, true);
            
            if (clusterSizes == null || clusterSizes.length == 0) {
                return -1;
            }
            
            writeIntArrayToFile(clusteringSizesFile, clusterSizes);
        }
        else {
            if (clusteringSizesFile.exists()) {
                clusterSizes = readIntArrayFromFile(clusteringSizesFile);
                
                if (clusterSizes == null || clusterSizes.length == 0) {
                    return -1;
                }
            }
            else {
                clusterSizes = mapWithMinimapAndExtractClusters(readsPath, targetPath, clustersdir,
                    numThreads, true, minimapOptions, stranded, maxEdgeClip,
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
                    minSeqDepth, usePacBioPreset, true);
                
                if (clusterSizes == null || clusterSizes.length == 0) {
                    return -1;
                }
                
                writeIntArrayToFile(clusteringSizesFile, clusterSizes);
            }
        }
                
        minSeqDepth = 1;
        cutRevCompArtifact = false;
        
        int numClusters = clusterSizes.length-1;
        System.out.println("Processing " + numClusters + " clusters...");
        int batchSize = 100;
        if (batchSize < numClusters * 0.1f) {
            batchSize = (int) Math.ceil(numClusters * 0.1f);
        }
        for (int cid=1; cid<=numClusters; ++cid) {
            int numReads = clusterSizes[cid];
            
            String clusterPrefix = clustersdir + File.separator + cid + File.separator + cid;
            
            File assemblyDoneStamp = new File(clusterPrefix + ".DONE");
            
            if (cid % batchSize == 0) {
                System.out.println(convertToRoundedPercent(cid/(float)numClusters) + "%: #" + cid);
            }
            
            if (forceOverwrite || !assemblyDoneStamp.exists()) {
                String inFastaPath = clusterPrefix + FASTA_EXT + GZIP_EXTENSION;
//                String subFastaPath = clusterPrefix + "_sub" + FASTA_EXT + GZIP_EXTENSION;
                String finalFastaPath = clusterPrefix + "_transcripts" + FASTA_EXT + GZIP_EXTENSION;
                String tmpPrefix = clusterPrefix + ".tmp";
                
                if (numReads == 1) {
                    symlinkRemoveExisting(inFastaPath, finalFastaPath);
                    touch(assemblyDoneStamp);
                }
                else {
//                    if (numReads >= minSeqDepth * 10) {
//                        // subsample if the cluster has too many reads
//                        ArrayList<BitSequence> bitSeqs = getListFromFile(inFastaPath);
//                        //Collections.sort(bitSeqs);
//                        kmerBased(bitSeqs, subFastaPath, bfSize, k, numHash, stranded,
//                                minSeqDepth, maxEdgeClip, false);
//                    }
//                    else {
//                        subFastaPath = inFastaPath;
//                    }

                    boolean status = uniqueOLC(inFastaPath, inFastaPath,
                                            finalFastaPath, tmpPrefix, null,
                                            numThreads, stranded, minimapOptions, maxEdgeClip,
                                            minAlnId, minOverlapMatches, maxIndelSize,
                                            minSeqDepth, usePacBioPreset, false, polyAReadNamesPath,
                                            polyaFinder, "rnabloom", false);

                    if (status) {
                        touch(assemblyDoneStamp);
                    }
                    else {
                        throw new IOException("Error processing cluster " + cid);
                    }
                }
            }
            else {
                System.out.println("WARNING: cluster " + cid + " was already processed!");
            }
        }
        
        return numClusters;
    }
    
    public static void main(String[] args) {
        //debug
        String seqFastaPath = "seed.fa";
        String uniqueFastaPath = "out.fa";
        int numThreads = 1;
        boolean align = true;
        String minimapOptions = "-K 250M -e 25 -f 0.0001 -c";
        boolean stranded = false;
        int maxEdgeClip = 50;
        float minAlnId = 0.6f;
        int minOverlapMatches = 150;
        int maxIndelSize = 50;
        int minSeqDepth = 2;
        boolean usePacBioPreset = false;
        boolean verbose = false;
        PolyATailFinder polyaFinder = new PolyATailFinder();
        polyaFinder.setProfile(PolyATailFinder.Profile.ONT);
        polyaFinder.setWindow(0);
        overlapWithMinimapAndExtractUnique(seqFastaPath, uniqueFastaPath,
            numThreads, align, minimapOptions, stranded,
            maxEdgeClip, minAlnId, minOverlapMatches, maxIndelSize,
            minSeqDepth, usePacBioPreset, verbose, polyaFinder);
    }
}
