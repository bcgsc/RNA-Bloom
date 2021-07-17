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
import static rnabloom.util.SeqUtils.filterFastq;
import static rnabloom.util.SeqUtils.longestSeq;
import static rnabloom.util.SeqUtils.reverseComplement;

/**
 *
 * @author Ka Ming Nip
 */
public class FastqFilteredSequenceIterator {
    private final Pattern seqPattern;
    private final Pattern qualPattern;
    private final int minAvgBaseQual;
    private final String[] fastqPaths;
    private final boolean reverseComplement;
    private int fileCursor = 0;
    private FastqReader reader = null;
    private final FastqRecord record = new FastqRecord();
    private boolean hasNext = true;
    
    public FastqFilteredSequenceIterator(String[] fastqPaths, Pattern seqPattern, Pattern qualPattern,
            int minAvgBaseQual, boolean reverseComplement) throws IOException {
        this.seqPattern = seqPattern;
        this.qualPattern = qualPattern;
        this.minAvgBaseQual = minAvgBaseQual;
        this.fastqPaths = fastqPaths;
        this.reverseComplement = reverseComplement;
        
        hasNext = fileCursor < fastqPaths.length;
        
        setReader(fastqPaths[fileCursor]);
    }
    
    private void setReader(String path) throws IOException {
        if (FastqReader.isCorrectFormat(path)) {
            reader = minAvgBaseQual > 0 ? new FastqFilteredReader(path, minAvgBaseQual) : new FastqReader(path);
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
            
            if (++fileCursor >= fastqPaths.length) {
                hasNext = false;
                return false;
            }
            
            hasNext = true;
            setReader(fastqPaths[fileCursor]);
            
            return this.hasNext();
        }
        
        return hasNext;
    }

    public synchronized String next() throws FileFormatException, IOException {
        if (!hasNext) {
            return null;
        }
        
        try {
            reader.nextWithoutName(record);
                        
            String seq = longestSeq(record.seq, record.qual, seqPattern, qualPattern);

            if (reverseComplement) {
                seq = reverseComplement(seq);
            }

            return seq;
        }
        catch (NoSuchElementException e) {
            reader.close();
            
            if (++fileCursor >= fastqPaths.length) {
                hasNext = false;
                return null;
            }
            
            hasNext = true;
            setReader(fastqPaths[fileCursor]);
            
            return this.next();
        }
    }
    
    public synchronized ArrayList<String> nextSegments() throws FileFormatException, IOException {
        if (!hasNext) {
            return null;
        }
        
        try {
            reader.nextWithoutName(record);
                        
            ArrayList<String> segments = filterFastq(record, seqPattern, qualPattern);

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
            
            if (++fileCursor >= fastqPaths.length) {
                hasNext = false;
                return null;
            }
            
            hasNext = true;
            setReader(fastqPaths[fileCursor]);
            
            return this.nextSegments();
        }
    }
}
