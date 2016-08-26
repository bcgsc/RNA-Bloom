/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.BufferedReader;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 *
 * @author kmnip
 */
public class FastaReader implements Iterator<String> {
    private final Iterator<String> itr;

    public FastaReader(BufferedReader br) {
        itr = br.lines().iterator();
    }

    @Override
    public boolean hasNext() {
        return itr.hasNext();
    }

    @Override
    public String next() {
        if (itr.hasNext()){
            if (! itr.next().startsWith(">")) {
                throw new NoSuchElementException("Line 1 of FASTA record is expected to start with '>'");
            }
            return itr.next();
        }
        
        throw new NoSuchElementException("End of file");
    }
}
