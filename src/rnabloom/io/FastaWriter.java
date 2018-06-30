/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package rnabloom.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Ka Ming Nip
 */
public class FastaWriter {
    private final BufferedWriter out;
    //private FileLock lock = null;
    
    public FastaWriter(String path, boolean append) throws IOException {
        out = new BufferedWriter(new FileWriter(path, append));
        /*
        FileOutputStream stream = new FileOutputStream(path, append);
        out = new BufferedWriter(new OutputStreamWriter(stream));
        
        FileChannel channel = stream.getChannel();
        
        for (lock = channel.tryLock(); lock == null;) {
            lock = channel.tryLock();
        }
        */
    }
    
    public synchronized void write(String header, String seq) throws IOException {
        out.write(">");
        out.write(header);
        out.newLine();
        out.write(seq);
        out.newLine();
    }
    
    public void close() throws IOException {
        //lock.release();
        out.flush();
        out.close();
    }
}
