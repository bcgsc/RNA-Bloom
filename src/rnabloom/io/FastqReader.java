/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
//import java.util.NoSuchElementException;
//import java.util.function.Supplier;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author kmnip
 */
public final class FastqReader {
    private final static String GZIP_EXTENSION = ".gz";
    private final static Pattern RECORD_NAME_PATTERN = Pattern.compile("^@([^\\s/]+)(?:/[12])?.*$");
    private final BufferedReader br;
    private final Iterator<String> itr;
    
    public FastqReader(String path) throws IOException {
        if (path.endsWith(GZIP_EXTENSION)) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(path)));
        }
        itr = br.lines().iterator();
    }

    public static boolean isFastq(String path) {
        try {
            // try to read the first line as a FASTQ file
            FastqReader reader = new FastqReader(path);
            reader.nextWithoutName(new FastqRecord());
        }
        catch (Exception e) {
            return false;
        }
        
        return true;
    }
    
    public boolean hasNext() {
        return itr.hasNext();
    }
    
    public void nextWithoutName(FastqRecord fr) throws Exception {
        String line1, line3;
        
        synchronized(this) {
            line1 = itr.next();
            fr.seq = itr.next();
            line3 = itr.next();
            fr.qual = itr.next();
        }
        
        if (line1.charAt(0) != '@') {
            throw new Exception("Line 1 of FASTQ record is expected to start with '@'");
        }

        if (line3.charAt(0) != '+') {
            throw new Exception("Line 3 of FASTQ record is expected to start with '+'");
        }
    }
    
    public void nextWithName(FastqRecord fr) throws Exception {
        String line1, line3;
        
        synchronized(this) {
            line1 = itr.next();
            fr.seq = itr.next();
            line3 = itr.next();
            fr.qual = itr.next();
        }
        
        Matcher m = RECORD_NAME_PATTERN.matcher(line1);
        if (m.matches()) {
            fr.name = m.group(1);
        }
        else {
            throw new Exception("Line 1 of a FASTQ record is expected to start with '@'");
        }

        if (line3.charAt(0) != '+') {
            throw new Exception("Line 3 of a FASTQ record is expected to start with '+'");
        }
    }
        
    public void close() throws IOException {
        br.close();
    }
}
