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
                // error
                return fr;
            }
            
            fr.seq = itr.next();
            
            String comment = itr.next();
            if (! comment.startsWith("+")) {
                // error
                return fr;
            }
            
            fr.qual = itr.next();
            
            return fr;
        }
        else {
            throw new NoSuchElementException();
        }
    }
}
