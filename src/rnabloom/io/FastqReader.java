/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.BufferedReader;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Supplier;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author kmnip
 */
public final class FastqReader implements Iterator<FastqRecord> {
    private final Iterator<String> itr;
    private final static Pattern RECORD_NAME_PATTERN = Pattern.compile("^@([^\\s/]+)(?:/[12])?.*$");
    private final Supplier<FastqRecord> nextFunction;

    public FastqReader(BufferedReader br, boolean storeReadName) {
        itr = br.lines().iterator();
        if (storeReadName) {
            nextFunction = this::nextWithName;
        }
        else {
            nextFunction = this::nextWithoutName;
        }
    }

    @Override
    public boolean hasNext() {
        return itr.hasNext();
    }

    @Override
    public FastqRecord next() {
        return nextFunction.get();
    }
    
    public FastqRecord nextWithoutName() {
        if (itr.hasNext()){
            FastqRecord fr = new FastqRecord();
            
            if (itr.next().charAt(0) != '@') {
                throw new NoSuchElementException("Line 1 of FASTQ record is expected to start with '@'");
            }
            
            fr.seq = itr.next();
            
            if (! itr.next().startsWith("+")) {
                throw new NoSuchElementException("Line 3 of FASTQ record is expected to start with '+'");
            }
            
            fr.qual = itr.next();
            
            return fr;
        }
        throw new NoSuchElementException("End of file");
    }
    
    public FastqRecord nextWithName() {
        if (itr.hasNext()){
            FastqRecord fr = new FastqRecord();
            
            Matcher m = RECORD_NAME_PATTERN.matcher(itr.next());
            if (m.matches()) {
                fr.name = m.group(1);
            }
            else {
                throw new NoSuchElementException("Line 1 of FASTQ record is expected to start with '@'");
            }
            
            fr.seq = itr.next();
            
            if (! itr.next().startsWith("+")) {
                throw new NoSuchElementException("Line 3 of FASTQ record is expected to start with '+'");
            }
            
            fr.qual = itr.next();
            
            return fr;
        }
        throw new NoSuchElementException("End of file");
    }
}
