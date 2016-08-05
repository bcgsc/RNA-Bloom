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
public class FastaReader implements Supplier<String> {
    private final Iterator<String> itr;

    public FastaReader(BufferedReader br) {
        itr = br.lines().iterator();
    }
    
    @Override
    public String get() {
        if (itr.hasNext()){
            if (! itr.next().startsWith(">")) {
                throw new NoSuchElementException("Line 1 of FASTA record is expected to start with '>'");
            }
            return itr.next();
        }
        else {
            throw new NoSuchElementException();
        }
    }
}
