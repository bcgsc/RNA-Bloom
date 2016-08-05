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

/**
 *
 * @author kmnip
 */
public final class FastqReader implements Supplier<FastqRecord> {
    private final Iterator<String> itr;

    public FastqReader(BufferedReader br) {
        itr = br.lines().iterator();
    }
    
    @Override
    public FastqRecord get() {
        if (itr.hasNext()){
            FastqRecord fr = new FastqRecord();
            
            String name = itr.next();
            if (! name.startsWith("@")) {
                throw new NoSuchElementException("Line 1 of FASTQ record is expected to start with '@'");
            }
            
            fr.seq = itr.next();
            
            if (! itr.next().startsWith("+")) {
                throw new NoSuchElementException("Line 3 of FASTQ record is expected to start with '+'");
            }
            
            fr.qual = itr.next();
            
            return fr;
        }
        else {
            throw new NoSuchElementException();
        }
    }
}
