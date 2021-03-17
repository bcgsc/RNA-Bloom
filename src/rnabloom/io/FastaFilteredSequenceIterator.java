/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.NoSuchElementException;
import java.util.regex.Pattern;
import static rnabloom.util.SeqUtils.filterFasta;
import static rnabloom.util.SeqUtils.longestSeq;
import static rnabloom.util.SeqUtils.reverseComplement;

/**
 *
 * @author Ka Ming Nip
 */
public class FastaFilteredSequenceIterator {
    private final Pattern seqPattern;
    private final String[] fastaPaths;
    private final boolean reverseComplement;
    private int fileCursor = 0;
    private FastaReader reader = null;
    private boolean hasNext = true;
    
    public FastaFilteredSequenceIterator(String[] fastaPaths, Pattern seqPattern, boolean reverseComplement) throws IOException {
        this.seqPattern = seqPattern;
        this.fastaPaths = fastaPaths;
        this.reverseComplement = reverseComplement;
        
        hasNext = fileCursor < fastaPaths.length;
        
        setReader(fastaPaths[fileCursor]);
    }
    
    private void setReader(String path) throws IOException {
        if (FastaReader.isCorrectFormat(path)) {
            reader = new FastaReader(path);
        }
        else {
            throw new FileFormatException("Incompatible file format for `" + path + "`");
        }
        
        System.out.println("Parsing `" + path + "`...");
    }
    
    public synchronized boolean hasNext() throws IOException {
        if (!hasNext) {
            return false;
        }

        hasNext = reader.hasNext();
        
        if (!hasNext) {
            reader.close();
            
            if (++fileCursor >= fastaPaths.length) {
                hasNext = false;
                return false;
            }
            
            hasNext = true;
            setReader(fastaPaths[fileCursor]);
            
            return this.hasNext();
        }
        
        return hasNext;
    }

    public synchronized String next() throws FileFormatException, IOException {
        if (!hasNext) {
            return null;
        }
        
        try {
            String seq = longestSeq(reader.next(), seqPattern);

            if (reverseComplement) {
                seq = reverseComplement(seq);
            }

            return seq;
        }
        catch (NoSuchElementException e) {
            reader.close();
            
            if (++fileCursor >= fastaPaths.length) {
                hasNext = false;
                return null;
            }
            
            hasNext = true;
            setReader(fastaPaths[fileCursor]);
            
            return this.next();
        }
    }
    
    public synchronized ArrayList<String> nextSegments() throws FileFormatException, IOException {
        if (!hasNext) {
            return null;
        }
        
        try {
            ArrayList<String> segments = filterFasta(reader.next(), seqPattern);

            if (reverseComplement) {
                int numRightSegments = segments.size();

                if (numRightSegments > 1) {
                    Collections.reverse(segments);
                }

                for (int i=0; i<numRightSegments; ++i) {
                    segments.set(i, reverseComplement(segments.get(i)));
                }
            }

            return segments;
        }
        catch (NoSuchElementException e) {
            reader.close();
            
            if (++fileCursor >= fastaPaths.length) {
                hasNext = false;
                return null;
            }
            
            hasNext = true;
            setReader(fastaPaths[fileCursor]);
            
            return this.nextSegments();
        }
    }
}
