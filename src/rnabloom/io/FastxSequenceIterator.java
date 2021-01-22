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

import java.io.IOException;
import java.util.NoSuchElementException;

/**
 *
 * @author Ka Ming Nip
 */
public class FastxSequenceIterator {
    private final String[] fastxPaths;
    private int fileCursor;
    private FastxReaderInterface reader;
    private String[] n = null;
    
    public FastxSequenceIterator(String[] fastxPaths) throws IOException {
        this.fastxPaths = fastxPaths;
        fileCursor = 0;
        setReader(fastxPaths[fileCursor]);
        if (hasNext()) {
            n = readNext();
        }
    }
    
    private void setReader(String path) throws IOException {
        if (FastqReader.isCorrectFormat(path)) {
            reader = new FastqReader(path);
        }
        else if (FastaReader.isCorrectFormat(path)) {
            reader = new FastaReader(path);
        }
        else {
            throw new FileFormatException("Incompatible file format for `" + path + "`");
        }
        
        System.out.println("Parsing `" + path + "`...");
    }
    
    public boolean hasNext() throws IOException {
        boolean hasNext = reader.hasNext();
        
        if (!hasNext) {
            reader.close();
            
            if (++fileCursor >= fastxPaths.length) {
                return false;
            }
            
            setReader(fastxPaths[fileCursor]);
            
            return this.hasNext();
        }
        
        return hasNext;
    }
    
    private String[] readNext() throws FileFormatException, IOException {        
        if (reader.hasNext()) {
            return reader.nextWithName();
        }
        else {
            reader.close();
            
            if (++fileCursor >= fastxPaths.length) {
                return null;
            }
            
            setReader(fastxPaths[fileCursor]);
            
            return this.readNext();
        }
    }
    
    public synchronized String[] next() throws IOException {
        if (n == null) {
            return null;
        }
        
        String[] p = this.n;
        this.n = readNext();
        
        return p;
    }
}
