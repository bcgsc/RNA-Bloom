/* 
 * Copyright (C) 2018-present BC Cancer Genome Sciences Centre
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
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import static rnabloom.io.Constants.BUFFER_SIZE;
import static rnabloom.io.Constants.GZIP_EXT;

/**
 *
 * @author Ka Ming Nip
 */
public class FastaWriter {
    private final BufferedWriter out;
    //private FileLock lock = null;
    
    public FastaWriter(String path, boolean append) throws IOException {
        if (path.toLowerCase().endsWith(GZIP_EXT)) {
            out = new BufferedWriter(
                    new OutputStreamWriter(
                        new FastGZIPOutputStream(
                            new FileOutputStream(path, append),
                            BUFFER_SIZE),
                        StandardCharsets.UTF_8),
                    BUFFER_SIZE);
        }
        else {
            out = new BufferedWriter(new FileWriter(path, append), BUFFER_SIZE);
        }
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
        out.write('>');
        out.write(header);
        out.write('\n');
        out.write(seq);
        out.write('\n');
    }
    
    public void close() throws IOException {
        //lock.release();
        out.flush();
        out.close();
    }
}
