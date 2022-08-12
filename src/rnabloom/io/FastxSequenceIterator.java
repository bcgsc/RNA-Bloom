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

import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;
import static rnabloom.util.SeqUtils.filterFasta;
import static rnabloom.util.SeqUtils.filterFastq;

/**
 *
 * @author Ka Ming Nip
 */
public final class FastxSequenceIterator {
    private final int minAvgBaseQual;
    private final String[] fastxPaths;
    private int fileCursor;
    private FastxReaderInterface reader;
    private String[] n = null;
    private final Pattern seqPattern;
    private final Pattern qualPattern;
    boolean isFastq = false;
    boolean isFasta = false;
    boolean removeNameSuffix = false;
    
    public FastxSequenceIterator(String[] fastxPaths, int minAvgBaseQual, boolean removeNameSuffix) throws IOException {
        this.fastxPaths = fastxPaths;
        this.minAvgBaseQual = minAvgBaseQual;
        fileCursor = 0;
        setReader(fastxPaths[fileCursor]);
        if (hasNext()) {
            n = readNext();
        }
        
        seqPattern = null;
        qualPattern = null;
        this.removeNameSuffix = removeNameSuffix;
    }
    
    public FastxSequenceIterator(String[] fastxPaths, int minAvgBaseQual, Pattern seqPattern, Pattern qualPattern) throws IOException {
        this.fastxPaths = fastxPaths;
        this.minAvgBaseQual = minAvgBaseQual;
        fileCursor = 0;
        setReader(fastxPaths[fileCursor]);
        if (hasNext()) {
            n = readNext();
        }
        
        this.seqPattern = seqPattern;
        this.qualPattern = qualPattern;
    }
    
    private void setReader(String path) throws IOException {
        if (FastqReader.isCorrectFormat(path)) {
            FastqReader fqReader = minAvgBaseQual > 0 ? new FastqFilteredReader(path, minAvgBaseQual) : new FastqReader(path);
            fqReader.setRemoveNameSuffix(removeNameSuffix);
            reader = fqReader;
            isFastq = true;
            isFasta = false;
        }
        else if (FastaReader.isCorrectFormat(path)) {
            FastaReader faReader = new FastaReader(path);
            faReader.setRemoveNameSuffix(removeNameSuffix);
            reader = faReader;
            isFasta = true;
            isFastq = false;
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
    
    public synchronized ArrayList<String> nextSegments() throws IOException {
        if (n == null) {
            return null;
        }
        
        String[] p = this.n;
        ArrayList<String> segments = null;
        if (isFasta) {
            segments = filterFasta(p[1], seqPattern);
        }
        else if (isFastq) {
            segments = filterFastq(p[1], p[2], seqPattern, qualPattern);
        }
        else {
            throw new FileFormatException("Incompatible file format");
        }
        
        this.n = readNext();
        
        return segments;
    }
}
