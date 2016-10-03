/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author kmnip
 */
public class FastaWriter {
    private final BufferedWriter out;
    
    public FastaWriter(String path, boolean append) throws IOException {
        out = new BufferedWriter(new FileWriter(path, append));
    }
    
    public void write(String header, String seq) throws IOException {
        out.write(">" + header);
        out.newLine();
        out.write(seq);
        out.newLine();
    }
    
    public void close() throws IOException {
        out.close();
    }
}