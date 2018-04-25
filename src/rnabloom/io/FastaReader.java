/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
//import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author kmnip
 */
public class FastaReader {
    private final static String GZIP_EXTENSION = ".gz";
    private final static Pattern RECORD_NAME_PATTERN = Pattern.compile("^>([^\\s/]+)(?:/[12])?.*$");
    private final Iterator<String> itr;
    private final BufferedReader br;
    
    public FastaReader(String path) throws IOException {
        if (path.endsWith(GZIP_EXTENSION) || Files.probeContentType(Paths.get(path)).equals("application/gzip")) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(path)));
        }
        itr = br.lines().iterator();
    }
    
    public FastaReader(File f) throws IOException {
        if (f.getName().endsWith(GZIP_EXTENSION)) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(f)));
        }
        itr = br.lines().iterator();
    }
    
    public static boolean isFasta(String path) {
        try {
            // try to get the first FASTA record
            FastaReader reader = new FastaReader(path);
            reader.next();
            reader.close();
        }
        catch (Exception e) {
            return false;
        }
        
        return true;
    }

//    @Override
    public boolean hasNext() {
        return itr.hasNext();
    }

//    @Override
    public synchronized String next() throws FileFormatException {
        if (itr.next().charAt(0) != '>') {
            throw new FileFormatException("Line 1 of a FASTA record is expected to start with '>'");
        }
        return itr.next();
    }

    public void nextWithName(FastaRecord fr) throws FileFormatException {
        String line1;
        
        synchronized(this) {
            line1 = itr.next();
            fr.seq = itr.next();
        }
    
        Matcher m = RECORD_NAME_PATTERN.matcher(line1);
        if (m.matches()) {
            fr.name = m.group(1);
        }
        else {
            throw new FileFormatException("Line 1 of a FASTA record is expected to start with '>'");
        }
    }
    
    public void close() throws IOException {
        br.close();
    }
}
