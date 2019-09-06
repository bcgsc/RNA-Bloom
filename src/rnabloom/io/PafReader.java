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
import java.util.zip.GZIPInputStream;

/**
 *
 * @author kmnip
 */
public class PafReader {
    private final static String GZIP_EXTENSION = ".gz";
    private final Iterator<String> itr;
    private final BufferedReader br;
    
    public PafReader(String path) throws IOException {
        if (path.toLowerCase().endsWith(GZIP_EXTENSION)) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))));
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(path)));
        }
        itr = br.lines().iterator();
    }
        
    public boolean hasNext() {
        return itr.hasNext();
    }
    
    public PafRecord next() {
        return new PafRecord(itr.next().trim().split("\t"));
    }
    
    public void next(PafRecord record) {
        String[] cols = itr.next().trim().split("\t");
        record.update(cols);
    }
    
    public void close() throws IOException {
        br.close();
    }
}
