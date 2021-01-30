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
import java.lang.ProcessBuilder.Redirect;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;
import static rnabloom.io.Constants.FASTA_EXT;
import rnabloom.io.FastaReader;

/**
 *
 * @author kmnip
 */
public class OverlapLayoutConcensus {
    
    private final static String GZIP_EXTENSION = ".gz";
    private final static String LOG_EXTENSION = ".log";
    
    private final static String MINIMAP2 = "minimap2";
    private final static String RACON = "racon";
    
    private final static String PRESET_PACBIO = "pb";
    private final static String PRESET_ONT = "ont";
    
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
    
    public static int[] overlapWithMinimapAndExtractClusters(String seqFastaPath, String clusterdir,
            long numReads, int numThreads, boolean align, String minimapOptions, boolean stranded,
            int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset) {
        
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        
        if (align) {
            minimapOptions = "-c " + minimapOptions;
        }
        
        String preset = usePacBioPreset ? PRESET_PACBIO : PRESET_ONT;
        command.add(MINIMAP2 + " -X -x ava-" + preset + " " + minimapOptions + " -t " + numThreads + " " + seqFastaPath + " " + seqFastaPath);
        
        int[] clusterSizes = null;
        try {            
            ProcessBuilder pb = new ProcessBuilder(command);

            File logFile = new File(clusterdir + LOG_EXTENSION);
            pb.redirectError(Redirect.to(logFile));
            
            Process process = pb.start();

            clusterSizes = extractClusters(seqFastaPath, numReads, process.getInputStream(), clusterdir, stranded, maxEdgeClip, 
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth);
            
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
    
    public static boolean overlapWithMinimapAndLayout(String seqFastaPath, String layoutFastaPath,
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

            if (!layout(seqFastaPath, process.getInputStream(), layoutFastaPath, stranded, maxEdgeClip, 
                    minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth)) {
                return false;
            }
            
            int exitStatus = process.waitFor();
            return exitStatus == 0;
        }
        catch (IOException | InterruptedException e) {
            return false;
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
    
    public static int[] extractClusters(String seqFastaPath, long numReads, InputStream overlapPafInputStream, String outdir,
            boolean stranded, int maxEdgeClip, float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean cutRevCompArtifact, int minSeqDepth) {
        try {
            Layout myLayout = new Layout(seqFastaPath, overlapPafInputStream, stranded, maxEdgeClip, minAlnId, 
                    minOverlapMatches, maxIndelSize, cutRevCompArtifact, minSeqDepth);
            return myLayout.extractClusters(outdir, numReads);
        } catch (Exception ex) {
            ex.printStackTrace();
            return null;
        }
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
        
    public static boolean concensusWithRacon(String queryFastaPath, String targetFastaPath, String mappingPafPath, String concensusFastaPath, int numThreads) {
        ArrayList<String> command = new ArrayList<>();
        command.add("/bin/sh");
        command.add("-c");
        command.add(RACON + " -t " + numThreads + " " + queryFastaPath + " " + mappingPafPath + " " + targetFastaPath + " > " + concensusFastaPath);
        //--no-trimming -u
        return runCommand(command, concensusFastaPath + LOG_EXTENSION);
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
    
    private static void symlinkRemoveExisting(String target, String link) throws IOException {
        Path targetPath = Paths.get(target);
        Path linkPath = Paths.get(link);
        
        Files.deleteIfExists(linkPath);
        
        Files.createSymbolicLink(linkPath, linkPath.getParent().relativize(targetPath));
    }
    
    public static boolean overlapLayout(String readsPath, String tmpPrefix, String layoutPath,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset) throws IOException {
        
        if (hasOnlyOneSequence(readsPath)) {
            symlinkRemoveExisting(readsPath, layoutPath);
            return true;
        }
        
        boolean status = overlapWithMinimapAndLayout(readsPath, layoutPath,
            numThreads, true, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
            minSeqDepth, usePacBioPreset);
        
        if (!status) {
            // PAF is empty
            symlinkRemoveExisting(readsPath, layoutPath);
        }
        
        return true;
    }

    public static boolean overlapLayoutConcensus(String readsPath, String tmpPrefix, String concensusPath, 
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset) throws IOException {
        String backbonesFa = tmpPrefix + "_backbones.fa";
        String mapPaf = tmpPrefix + "_map.paf.gz";
        
        Files.deleteIfExists(FileSystems.getDefault().getPath(backbonesFa));
        Files.deleteIfExists(FileSystems.getDefault().getPath(mapPaf));
        
        if (hasOnlyOneSequence(readsPath)) {
            symlinkRemoveExisting(readsPath, concensusPath);
            return true;
        }
        
        boolean status = overlapWithMinimapAndLayout(readsPath, backbonesFa,
            numThreads, false, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
            minSeqDepth, usePacBioPreset);
        
        if (!status) {
            // PAF is empty
            symlinkRemoveExisting(readsPath, concensusPath);
            return true;
        }
        
        if (!mapWithMinimap(readsPath, backbonesFa, mapPaf, numThreads, minimapOptions, usePacBioPreset)) {
            return false;
        }
        
        return concensusWithRacon(readsPath, backbonesFa, mapPaf, concensusPath, numThreads);
    }
    
    public static boolean overlapLayoutConcensus2(String readsPath, String tmpPrefix, String concensusPath, 
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize,
            boolean usePacBioPreset, boolean alignReads, boolean alignBackbones) throws IOException {
        String backbonesFa = tmpPrefix + "_backbones.fa";
        //String polishedFa = tmpPrefix + "_polished.fa";
        String backbonesFa2 = tmpPrefix + "_backbones2.fa";
        String mapPaf = tmpPrefix + "_map.paf.gz";

        Files.deleteIfExists(FileSystems.getDefault().getPath(backbonesFa));
        //Files.deleteIfExists(FileSystems.getDefault().getPath(polishedFa));
        Files.deleteIfExists(FileSystems.getDefault().getPath(backbonesFa2));
        Files.deleteIfExists(FileSystems.getDefault().getPath(mapPaf));

        if (hasOnlyOneSequence(readsPath)) {
            symlinkRemoveExisting(readsPath, concensusPath);
            return true;
        }
        
        boolean status = overlapWithMinimapAndLayout(readsPath, backbonesFa,
            numThreads, alignReads, minimapOptions, stranded, maxEdgeClip,
            minAlnId, minOverlapMatches, maxIndelSize, false,
            1, usePacBioPreset);

        if (!status) {
            // either PAF is empty or no backbones can be made
            symlinkRemoveExisting(readsPath, concensusPath);
            return true;
        }

        if (hasOnlyOneSequence(backbonesFa)) {
            // polish backbone #1
            
            status = mapWithMinimap(readsPath, backbonesFa, mapPaf, numThreads, minimapOptions, usePacBioPreset);
            if (!status) {
                return false;
            }

            status = concensusWithRacon(readsPath, backbonesFa, mapPaf, concensusPath, numThreads);
        }
        else {
            // layout backbone #2

            status = overlapWithMinimapAndLayout(backbonesFa, backbonesFa2,
                        numThreads, alignBackbones, minimapOptions, stranded, maxEdgeClip,
                        minAlnId, minOverlapMatches, maxIndelSize, false,
                        1, usePacBioPreset);

            if (!status) {
                // either PAF is empty or no backbones can be made
//                symlinkRemoveExisting(backbonesFa2, concensusPath);
//                return true;
                backbonesFa2 = backbonesFa;
            }

            // polish backbone #2

            status = mapWithMinimap(readsPath, backbonesFa2, mapPaf, numThreads, minimapOptions, usePacBioPreset);
            if (!status) {
                return false;
            }

            status = concensusWithRacon(readsPath, backbonesFa2, mapPaf, concensusPath, numThreads);
        }
        
        return status;
    }
    
    public static int clusteredOLC(String readsPath, String clustersdir, long numReads,
            int numThreads, boolean stranded, String minimapOptions, int maxEdgeClip,
            float minAlnId, int minOverlapMatches, int maxIndelSize, boolean cutRevCompArtifact,
            int minSeqDepth, boolean usePacBioPreset) throws IOException {
        
        if (!Files.isDirectory(FileSystems.getDefault().getPath(clustersdir))) {
            Files.createDirectory(FileSystems.getDefault().getPath(clustersdir));
        }
                
        int[] clusterSizes = overlapWithMinimapAndExtractClusters(readsPath, clustersdir,
                            numReads,
                            numThreads, false, minimapOptions, stranded, maxEdgeClip,
                            minAlnId, minOverlapMatches, maxIndelSize, cutRevCompArtifact,
                            minSeqDepth, usePacBioPreset);
        
        if (clusterSizes == null) {
            return -1;
        }
        
        int numClusters = clusterSizes.length;
        
        for (int c=0; c<numClusters; ++c) {
            // do not align during AVA overlap if the cluster has to many reads
            boolean alignReads = clusterSizes[c] <= 100000;
            boolean alignBackbones = alignReads;
            
            int cid = c+1;
            System.out.println("Processing cluster #" + cid + " of " + numClusters + "...");
            
            String clusterPrefix = clustersdir + File.separator + cid;
            String clusterPath = clusterPrefix + FASTA_EXT;
            String tmpPrefix = clusterPrefix + "_tmp";
            String concensusPath = clusterPrefix + "_transcripts" + FASTA_EXT;
            
            boolean status = overlapLayoutConcensus2(clusterPath, tmpPrefix, concensusPath, 
                            numThreads, stranded, minimapOptions, maxEdgeClip,
                            minAlnId, minOverlapMatches, maxIndelSize,
                            usePacBioPreset, alignReads, alignBackbones);
            if (!status) {
                throw new IOException("Error processing cluster " + cid);
            }
        }
        
        return numClusters;
    }
    
    public static void main(String[] args) {
        //debug
    }
}
