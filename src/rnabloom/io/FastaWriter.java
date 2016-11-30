/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.channels.FileChannel;
import java.nio.channels.FileLock;

/**
 *
 * @author kmnip
 */
public class FastaWriter {
    //private final BufferedWriter out;
    private OutputStreamWriter out;
    private FileLock lock = null;
    
    public FastaWriter(String path, boolean append) throws IOException {
        //out = new BufferedWriter(new FileWriter(path, append));
        FileOutputStream stream = new FileOutputStream(path, append);
        out = new OutputStreamWriter(stream);
        
        FileChannel channel = stream.getChannel();
        
        while (lock == null) {
            lock = channel.tryLock();
        }
    }
    
    public void write(String header, String seq) throws IOException {
        out.write(">");
        out.write(header);
        out.write("\n");
        out.write(seq);
        out.write("\n");
    }
    
    public void close() throws IOException {
        lock.release();
        out.close();
    }
}
