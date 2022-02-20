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
package rnabloom.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.zip.GZIPInputStream;
import static rnabloom.io.Constants.BUFFER_SIZE;
import static rnabloom.io.Constants.GZIP_EXT;
import rnabloom.io.FastGZIPOutputStream;
import rnabloom.io.FastaReader;

/**
 *
 * @author Ka Ming Nip
 */
public class FileUtils {
    
    public static BufferedReader getTextFileReader(String path) throws FileNotFoundException, IOException {
        if (path.toLowerCase().endsWith(GZIP_EXT)) {
            return new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path), BUFFER_SIZE)), BUFFER_SIZE);
        }
        else {
            return new BufferedReader(new InputStreamReader(new FileInputStream(path)), BUFFER_SIZE);
        }
    }
    
    public static Writer getTextFileWriter(String path, boolean append) throws FileNotFoundException, IOException {
        if (path.toLowerCase().endsWith(GZIP_EXT)) {
            return new OutputStreamWriter(
                        new FastGZIPOutputStream(
                            new FileOutputStream(path, append),
                            BUFFER_SIZE),
                        StandardCharsets.UTF_8);
        }
        else {
            return new BufferedWriter(new FileWriter(path, append), BUFFER_SIZE);
        }
    }
    
    public static void touch(File f) throws IOException {
        f.getParentFile().mkdirs();
        if (!f.createNewFile()){
            f.setLastModified(System.currentTimeMillis());
        }
    }
    
    public static void symlinkRemoveExisting(String target, String link) throws IOException {
        Path targetPath = Paths.get(target);
        Path linkPath = Paths.get(link);
        
        Files.deleteIfExists(linkPath);
        
        Files.createSymbolicLink(linkPath, linkPath.getParent().relativize(targetPath));
    }
    
    public static void deleteIfExists(String p) throws IOException {
        Files.deleteIfExists(FileSystems.getDefault().getPath(p));
    }
    
    public static void saveStringToFile(String path, String str) throws IOException {
        FileWriter writer = new FileWriter(path, false);
        writer.write(str);
        writer.close();
    }
    
    public static String loadStringFromFile(String path) throws FileNotFoundException, IOException {
        BufferedReader br = getTextFileReader(path);
        String line = br.readLine().strip();
        return line;
    }
    
    public static void writeIntArrayToFile(File f, int[] arr) throws IOException {
        BufferedWriter fr = new BufferedWriter(new FileWriter(f, false));
        for (int i : arr) {
            fr.write(Integer.toString(i));
            fr.write('\n');
        }
        fr.close();
    }
    
    public static int[] readIntArrayFromFile(File f) throws FileNotFoundException, IOException {
        ArrayDeque<Integer> tmpQueue = new ArrayDeque<>();
        BufferedReader br = new BufferedReader(new FileReader(f));
        for (String line; (line = br.readLine()) != null; ) {
            tmpQueue.add(Integer.parseInt(line.trim()));
        }
        
        br.close();
        
        return tmpQueue.stream().mapToInt(i->i).toArray();
    }
    
    public static boolean hasOnlyOneSequence(String fasta) throws IOException {
        FastaReader reader = new FastaReader(fasta);
        int numSeq = 0;
        while (reader.hasNext()) {
            reader.next();
            if (++numSeq > 1) {
                reader.close();
                return false;
            }
        }
        reader.close();
        
        return numSeq == 1;
    }

    public static void fileToStringCollection(String path, Collection<String> c) throws FileNotFoundException, IOException {
        BufferedReader br = getTextFileReader(path);
        
        for (String line; (line = br.readLine()) != null; ) {
            c.add(line.trim());
        }
        
        br.close();
    }
}
