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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import static rnabloom.io.Constants.FASTA_EXT;
import rnabloom.io.ExtendedPafRecord;
import static rnabloom.util.FileUtils.deleteIfExists;
import static rnabloom.util.FileUtils.hasOnlyOneSequence;
import static rnabloom.util.FileUtils.readIntArrayFromFile;
import static rnabloom.util.FileUtils.symlinkRemoveExisting;
import static rnabloom.util.FileUtils.touch;
import static rnabloom.util.FileUtils.writeIntArrayToFile;
import static rnabloom.util.PafUtils.hasGoodAlignment;
import static rnabloom.util.PafUtils.hasLargeOverlap;

/**
 *
 * @author kmnip
 */
public class OverlapLayoutConsensus {
    
    private final static String GZIP_EXTENSION = ".gz";
    private final static String LOG_EXTENSION = ".log";
    
    private final static String MINIMAP2 = "minimap2";
    private final static String RACON = "racon";
    
    private final static String PRESET_PACBIO = "pb";
    private final static String PRESET_ONT = "ont";
    
    public static enum STATUS {SUCCESS, EMPTY, FAIL};
    
    private static boolean runCommand(List<String> command, String logPath) {
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);
            if (logPath != null) {
                // write stdout and stderr to file to avoid hanging due to buffer overflow
                pb.redirectErrorStream(true);
                File logFile = new File(logPath);
                pb.redirectOutput(Redirect.to(logFile));
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
                pb.redirectError(Redirect.to(logFile));
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
            boolean align, String options, boolean usePacBioPreset) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align) {
            options = "-c " + options;
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -X -x ava-" + preset + " " + options + " -t " + numThreads + " " + seqFastaPath + " " + seqFastaPath + " | gzip -c > " + outPafPath);
        
        return runCommand(command, outPafPath + LOG_EXTENSION);
    }
    
    public static STATUS overlapWithMinimapAndExtractUnique(String seqFastaPath, String uniqueFastaPath,
            int numThreads, boolean align, String minimapOptions, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            int minSeqDepth, boolean usePacBioPreset) {
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align) {
            minimapOptions = "-c " + minimapOptions;
        }
                
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -X -x ava-" + preset + " " + minimapOptions + " -t " + numThreads + " " + seqFastaPath + " " + seqFastaPath);
        
        STATUS status = STATUS.FAIL;
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(uniqueFastaPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            Layout myLayout = new Layout(seqFastaPath, process.getInputStream(), stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, false, minSeqDepth);
            long numUnique = myLayout.extractUnique(uniqueFastaPath);
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
    
    public static int[] overlapWithMinimapAndExtractClusters(String seqFastaPath, String clusterdir,
            int numThreads, boolean align, String minimapOptions, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, int maxMergedClusterSize) {
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align) {
            minimapOptions = "-c " + minimapOptions;
        }
        
        if (!minimapOptions.contains("-g ")) {
            minimapOptions += " -g " + 2 * maxIndelSize;
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -X -x ava-" + preset + " " + minimapOptions + " -t " + numThreads + " " + seqFastaPath + " " + seqFastaPath);
        
        int[] clusterSizes = null;
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(clusterdir + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();
            
            Layout myLayout = new Layout(seqFastaPath, process.getInputStream(), stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth);
            clusterSizes = myLayout.extractClusters(clusterdir, maxMergedClusterSize);
            
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
            int minSeqDepth, boolean usePacBioPreset) {
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align) {
            minimapOptions = "-c " + minimapOptions;
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -X -x ava-" + preset + " " + minimapOptions + " -t " + numThreads + " " + seqFastaPath + " " + seqFastaPath);
        
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(layoutFastaPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            boolean reduced = layout(seqFastaPath, process.getInputStream(), layoutFastaPath, stranded, maxEdgeClip, 
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth);
            
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
            int numThreads, String options, boolean usePacBioPreset) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");

        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x map-" + preset + " -c " + options + " -t " + numThreads + " " + targetFastaPath + " " + queryFastaPath + " | gzip -c > " + outPafPath);
        
        return runCommand(command, outPafPath + LOG_EXTENSION);
    }
    
    public static STATUS mapWithMinimapFiltered(String queryFastaPath, String targetFastaPath, String outPafPath,
            int numThreads, String options, boolean usePacBioPreset, boolean stranded, int maxIndelSize, int minOverlapMatches, float minAlnId) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");

        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -x map-" + preset + " -c " + options + " -t " + numThreads + " " + targetFastaPath + " " + queryFastaPath);

        try {
            Writer bw = new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outPafPath, false)), "UTF-8");
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(outPafPath + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            long numGoodRecords = 0;
            ExtendedPafRecord record = new ExtendedPafRecord();
            BufferedReader br = new BufferedReader(new InputStreamReader(process.getInputStream()));
            for (String line; (line = br.readLine()) != null; ) {
                record.update(line.trim().split("\t"));
                if (!stranded || !record.reverseComplemented) {
                    if (hasLargeOverlap(record, minOverlapMatches) &&
                            hasGoodAlignment(record, maxIndelSize, minAlnId)) {
                        bw.write(line);
                        bw.write('\n');
                        ++numGoodRecords;
                    }
                }
            }
            br.close();
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
        
        return STATUS.SUCCESS;
    }
        
    public static boolean layout(String seqFastaPath, InputStream overlapPafInputStream, String backboneFastaPath,
            boolean stranded, int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth) {
        try {
            Layout myLayout = new Layout(seqFastaPath, overlapPafInputStream, stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth);
            return myLayout.writeBackboneSequences(backboneFastaPath);
        } catch (Exception ex) {
            ex.printStackTrace();
            return false;
        }
    }
        
    public static boolean consensusWithRacon(String queryFastaPath, String targetFastaPath, 
            String mappingPafPath, String consensusFastaPath, int numThreads, boolean keepUnpolished) {
        String extraOptions = "";
        if (keepUnpolished) {
            extraOptions += " --no-trimming -u ";
        }
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        command.add(RACON + extraOptions + " -t " + numThreads + " " + queryFastaPath + " " + mappingPafPath + " " + targetFastaPath + " > " + consensusFastaPath);
        //--no-trimming -u
        return runCommand(command, consensusFastaPath + LOG_EXTENSION);
    }
        
    public static boolean overlapLayout(String readsPath, String layoutPath,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset) throws IOException {
        
        if (hasOnlyOneSequence(readsPath)) {
            symlinkRemoveExisting(readsPath, layoutPath);
            return true;
        }
        
        STATUS s = overlapWithMinimapAndLayout(readsPath, layoutPath,
            numThreads, true, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
            minSeqDepth, usePacBioPreset);
        
        switch (s) {
            case FAIL:
                return false;
            case EMPTY:
                symlinkRemoveExisting(readsPath, layoutPath);
        }
        
        return true;
    }

    public static boolean mapAndConsensus(String readsPath, String targetPath, String tmpPrefix, String consensusPath, 
            int numThreads, String minimapOptions, boolean usePacBioPreset, boolean keepUnpolished) throws IOException {
        String mapPaf = tmpPrefix + "_map.paf.gz";
        deleteIfExists(mapPaf);
        
        if (!mapWithMinimap(readsPath, targetPath, mapPaf, numThreads, minimapOptions, usePacBioPreset)) {
            return false;
        }
        
        return consensusWithRacon(readsPath, targetPath, mapPaf, consensusPath, numThreads, keepUnpolished);
    }
    
    public static boolean overlapLayoutConsensus(String readsPath, String tmpPrefix, String consensusPath, 
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, boolean align, boolean keepUnpolished) throws IOException {
        String backbonesFa = tmpPrefix + "_backbones.fa";
        String mapPaf = tmpPrefix + "_map.paf.gz";
        
        deleteIfExists(backbonesFa);
        deleteIfExists(mapPaf);
        
        STATUS s = overlapWithMinimapAndLayout(readsPath, backbonesFa,
            numThreads, align, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
            minSeqDepth, usePacBioPreset);
        
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
        
        return consensusWithRacon(readsPath, backbonesFa, mapPaf, consensusPath, numThreads, keepUnpolished);
    }
    
    public static boolean overlapLayoutConsensus2(String readsPath, String tmpPrefix, String consensusPath, 
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean usePacBioPreset, boolean alignReads, boolean alignBackbones, boolean keepUnpolished) throws IOException {
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
            1, usePacBioPreset);

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

            return consensusWithRacon(readsPath, backbonesFa, mapPaf, consensusPath, numThreads, keepUnpolished);
        }
        else {
            // layout backbone #2

            s = overlapWithMinimapAndLayout(backbonesFa, backbonesFa2,
                        numThreads, alignBackbones, minimapOptions, stranded, maxEdgeClip,
                        minAlnId, minOverlapMatches, maxIndelSize, false,
                        1, usePacBioPreset);

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

            return consensusWithRacon(readsPath, backbonesFa2, mapPaf, consensusPath, numThreads, keepUnpolished);
        }
    }
    
    public static boolean uniqueOLC(String readsPath, String tmpPrefix, String outFastaPath,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize,
            int minSeqDepth, boolean usePacBioPreset) throws IOException {

        String uniqueFastaPath = tmpPrefix + "1.nr" + FASTA_EXT + GZIP_EXTENSION;
        String readsToUniquePafPath = tmpPrefix + "2.map.paf.gz";
        String polishedUniqueFastaPath = tmpPrefix + "3.pol" + FASTA_EXT;
        String backboneFastaPath = tmpPrefix + "4.backbone" + FASTA_EXT + GZIP_EXTENSION;
        String readsToBackbonePafPath = tmpPrefix + "5.map.paf.gz";
        
        deleteIfExists(uniqueFastaPath);
        deleteIfExists(readsToUniquePafPath);
        deleteIfExists(polishedUniqueFastaPath);
        deleteIfExists(backboneFastaPath);
        deleteIfExists(readsToBackbonePafPath);
        deleteIfExists(outFastaPath);
        
//        String minimapOptionsNoGaps = minimapOptions;
//        if (!minimapOptionsNoGaps.contains("-g ")) {
//            minimapOptionsNoGaps += " -g " + 2 * maxIndelSize;
//        }
        
        // overlap all reads and extract unique reads
        STATUS status = overlapWithMinimapAndExtractUnique(readsPath, uniqueFastaPath,
            numThreads, false, minimapOptions, stranded,
            maxEdgeClip, minAlnId, minOverlapMatches, maxIndelSize,
            minSeqDepth, usePacBioPreset);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        // map all reads to unique reads
        status = mapWithMinimapFiltered(readsPath, uniqueFastaPath, readsToUniquePafPath,
            numThreads, minimapOptions, usePacBioPreset, stranded, maxIndelSize,
            minOverlapMatches, minAlnId);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        // derive concensus for unique reads
        boolean success = consensusWithRacon(readsPath, uniqueFastaPath, 
            readsToUniquePafPath, polishedUniqueFastaPath, numThreads, true);
        if (!success) {
            return false;
        }
        
        // overlap concensus unique reads and layout
        success = overlapLayout(polishedUniqueFastaPath, backboneFastaPath,
            numThreads, stranded, minimapOptions, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, false,
            minSeqDepth, usePacBioPreset);
        if (!success) {
            return false;
        }
        
        // map all reads to backbone
        status = mapWithMinimapFiltered(readsPath, backboneFastaPath, readsToBackbonePafPath,
            numThreads, minimapOptions, usePacBioPreset, stranded, maxIndelSize,
            minOverlapMatches, minAlnId);
        if (status != STATUS.SUCCESS) {
            return false;
        }
        
        // derive concensus for unique reads
        success = consensusWithRacon(readsPath, backboneFastaPath, 
            readsToBackbonePafPath, outFastaPath, numThreads, true);
        
        return success;
    }
    
    public static int clusteredOLC(String readsPath, String clustersdir,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset, int maxMergedClusterSize, boolean forceOverwrite) throws IOException {
        
        File clusteringSizesFile = new File(clustersdir + File.separator + "cluster_sizes.txt");
        int[] clusterSizes = null;
        
        if (forceOverwrite) {
            Files.deleteIfExists(clusteringSizesFile.toPath());
            
            clusterSizes = overlapWithMinimapAndExtractClusters(readsPath, clustersdir,
                    numThreads, false, minimapOptions, stranded, maxEdgeClip,
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
                    minSeqDepth, usePacBioPreset, maxMergedClusterSize);
            
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
                    minSeqDepth, usePacBioPreset, maxMergedClusterSize);
                
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
            
            File assemblyDoneStamp = new File(clustersdir + File.separator + cid + ".DONE");
            System.out.println("Processing cluster #" + cid + " of " + numClusters + "...");
            
            if (forceOverwrite || !assemblyDoneStamp.exists()) {
                String clusterPrefix = clustersdir + File.separator + cid;
                String inFastaPath = clusterPrefix + FASTA_EXT;
//                String tmpPrefix = clusterPrefix + "_tmp";
                String finalFastaPath = clusterPrefix + "_transcripts" + FASTA_EXT;

                // do not align during AVA overlap if the cluster has to many reads
                int numReads = clusterSizes[c];
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
                                    usePacBioPreset);
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
    
    public static void main(String[] args) {
        //debug
    }
}
