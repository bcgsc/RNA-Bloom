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
//    private final Supplier<FastqRecord> nextFunction;
    public final FastqRecord record = new FastqRecord();
    
    public FastqReader(String path) throws IOException {
        if (path.endsWith(GZIP_EXTENSION)) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(path)));
        }
        itr = br.lines().iterator();
//        if (storeReadName) {
//            nextFunction = this::nextWithName;
//        }
//        else {
//            nextFunction = this::nextWithoutName;
//        }
    }

//    @Override
    public boolean hasNext() {
        return itr.hasNext();
    }

//    @Override
//    public synchronized FastqRecord next() {
//        return nextFunction.get();
//    }
    
    public synchronized void nextWithoutNameFunction() throws Exception {
//        if (itr.hasNext()){
            if (itr.next().charAt(0) != '@') {
                throw new Exception("Line 1 of FASTQ record is expected to start with '@'");
            }
            
            record.seq = itr.next();
            
            if (itr.next().charAt(0) != '+') {
                throw new Exception("Line 3 of FASTQ record is expected to start with '+'");
            }
            
            record.qual = itr.next();
//        }
//        else {
//            throw new NoSuchElementException("End of file");
//        }
    }
    
    public synchronized FastqRecord nextWithoutName() throws Exception {
//        if (itr.hasNext()){
            FastqRecord fr = new FastqRecord();
            
            if (itr.next().charAt(0) != '@') {
                throw new Exception("Line 1 of FASTQ record is expected to start with '@'");
            }
            
            fr.seq = itr.next();
            
            if (itr.next().charAt(0) != '+') {
                throw new Exception("Line 3 of FASTQ record is expected to start with '+'");
            }
            
            fr.qual = itr.next();
            
            return fr;
//        }
//        throw new NoSuchElementException("End of file");
    }

    public synchronized void nextWithNameFunction() throws Exception {
//        if (itr.hasNext()){
            Matcher m = RECORD_NAME_PATTERN.matcher(itr.next());
            if (m.matches()) {
                record.name = m.group(1);
            }
            else {
                throw new Exception("Line 1 of a FASTQ record is expected to start with '@'");
            }
            
            record.seq = itr.next();
            
            if (itr.next().charAt(0) != '+') {
                throw new Exception("Line 3 of a FASTQ record is expected to start with '+'");
            }
            
            record.qual = itr.next();
//        }
//        else {
//            throw new NoSuchElementException("Reached the end of file");
//        }
    }
    
    public synchronized FastqRecord nextWithName() throws Exception {
//        if (itr.hasNext()){
            FastqRecord fr = new FastqRecord();
            
            Matcher m = RECORD_NAME_PATTERN.matcher(itr.next());
            if (m.matches()) {
                fr.name = m.group(1);
            }
            else {
                throw new Exception("Line 1 of a FASTQ record is expected to start with '@'");
            }
            
            fr.seq = itr.next();
            
            if (itr.next().charAt(0) != '+') {
                throw new Exception("Line 3 of a FASTQ record is expected to start with '+'");
            }
            
            fr.qual = itr.next();
            
            return fr;
//        }
//        throw new NoSuchElementException("Reached the end of file");
    }
        
    public void close() throws IOException {
        br.close();
    }
}
