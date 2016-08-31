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
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

/**
 *
 * @author kmnip
 */
public class FastaReader implements Iterator<String> {
    private final static String GZIP_EXTENSION = ".gz";
    private final Iterator<String> itr;
    private final BufferedReader br;
    
    public FastaReader(String path) throws IOException {
        if (path.endsWith(GZIP_EXTENSION)) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(path)));
        }
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
                throw new NoSuchElementException("Line 1 of a FASTA record is expected to start with '>'");
            }
            return itr.next();
        }
        
        throw new NoSuchElementException("Reached the end of file");
    }
    
    public void close() throws IOException {
        br.close();
    }
}
