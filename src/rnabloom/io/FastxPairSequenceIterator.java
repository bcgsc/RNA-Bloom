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
import java.util.regex.Pattern;

/**
 *
 * @author gengar
 */
public class FastxPairSequenceIterator {
    private final FastxFilePair[] fastxPairs;
    private int fileCursor;
    private final Pattern seqPattern;
    private final Pattern qualPattern;
    private FastxPairReader reader;
    private PairedReadSegments n = null;
    
    public FastxPairSequenceIterator(FastxFilePair[] fastxPairs, Pattern seqPattern, Pattern qualPattern) throws IOException {
        this.seqPattern = seqPattern;
        this.qualPattern = qualPattern;
        this.fastxPairs = fastxPairs;
        this.fileCursor = 0;
        setReader(fastxPairs[fileCursor]);
        if (hasNext()) {
            n = readNext();
        }
    }
    
    private void setReader(FastxFilePair fxPair) throws IOException {
        if (FastqReader.isCorrectFormat(fxPair.leftPath) && FastqReader.isCorrectFormat(fxPair.rightPath)) {
            reader = new FastqPairReader(fxPair.leftPath, fxPair.rightPath, qualPattern, seqPattern, fxPair.leftRevComp, fxPair.rightRevComp);
        }
        else if (FastaReader.isCorrectFormat(fxPair.leftPath) && FastaReader.isCorrectFormat(fxPair.rightPath)) {
            reader = new FastaPairReader(fxPair.leftPath, fxPair.rightPath, seqPattern, fxPair.leftRevComp, fxPair.rightRevComp);
        }
        else {
            throw new FileFormatException("Incompatible file format for `" + fxPair.leftPath + "` and `" + fxPair.rightPath + "`");
        }
        
        System.out.println("Parsing `" + fxPair.leftPath + "` and `" + fxPair.rightPath + "`...");
    }
    
    public boolean hasNext() throws IOException {
        boolean hasNext = reader.hasNext();
        
        if (!hasNext) {
            reader.close();
            
            if (++fileCursor >= fastxPairs.length) {
                return false;
            }
            
            setReader(fastxPairs[fileCursor]);
            
            return this.hasNext();
        }
        
        return hasNext;
    }
    
    private PairedReadSegments readNext() throws FileFormatException, IOException {        
        try {
            return reader.next();
        }
        catch (NoSuchElementException e) {
            reader.close();
            
            if (++fileCursor >= fastxPairs.length) {
                return null;
            }
            
            setReader(fastxPairs[fileCursor]);
            
            return this.readNext();
        }
    }
    
    public synchronized PairedReadSegments next() throws IOException {
        if (n == null) {
            return null;
        }
        
        PairedReadSegments p = this.n;
        this.n = readNext();
        
        return p;
    }
}
