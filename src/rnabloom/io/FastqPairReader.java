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
import java.util.ArrayList;
import java.util.Collections;
import java.util.regex.Pattern;
import static rnabloom.util.SeqUtils.filterFastq;
import static rnabloom.util.SeqUtils.reverseComplement;

/**
 *
 * @author Ka Ming Nip
 */
public class FastqPairReader implements FastxPairReader {
    private final FastqReader leftReader;
    private final FastqReader rightReader;
    private final Pattern qualPattern;
    private final Pattern seqPattern;
    private final NextFunction<PairedReadSegments> nextFunction;
                
    public FastqPairReader(String leftPath, String rightPath, Pattern qualPattern, Pattern seqPattern, boolean revCompLeft, boolean revCompRight) throws IOException {
        leftReader = new FastqReader(leftPath);
        rightReader = new FastqReader(rightPath);
        
        this.qualPattern = qualPattern;
        this.seqPattern = seqPattern;
        
        if (revCompLeft) { 
            if (revCompRight) {
                nextFunction = this::nextRR;
            }
            else {
                nextFunction = this::nextRF;
            }
        }
        else {
            if (revCompRight) {
                nextFunction = this::nextFR;
            }
            else {
                nextFunction = this::nextFF;
            }                
        }
    }

    public interface NextFunction<PairedReadSegments> {
        PairedReadSegments get() throws FileFormatException;
    }
    
    @Override
    public boolean hasNext() {
        return this.leftReader.hasNext() && this.rightReader.hasNext();
    }

    @Override
    public synchronized PairedReadSegments next() throws FileFormatException {
        return nextFunction.get();
    }

    private synchronized PairedReadSegments nextFF() throws FileFormatException {
        PairedReadSegments p = new PairedReadSegments();

        FastqRecord frLeft = new FastqRecord();
        FastqRecord frRight = new FastqRecord();
        
        try {
            leftReader.nextWithName(frLeft);
            rightReader.nextWithName(frRight);
        }
        catch (FileFormatException e) {
            throw new FileFormatException(e.getMessage());
        }
        
        if (!frLeft.name.equals(frRight.name)) {
            throw new FileFormatException("Inconsistent record names: \"" + frLeft.name + "\" and \"" + frRight.name + "\"");
        }

        p.left = filterFastq(frLeft, seqPattern, qualPattern);
        
//        int numBasesTrimmed = frLeft.seq.length();
//        p.originalLeftLength = numBasesTrimmed;
//        for (String s : p.left) {
//            numBasesTrimmed -= s.length();
//        }
//        p.numLeftBasesTrimmed = numBasesTrimmed;
        
        p.right = filterFastq(frRight, seqPattern, qualPattern);
        
//        numBasesTrimmed = frRight.seq.length();
//        p.originalRightLength = numBasesTrimmed;
//        for (String s : p.right) {
//            numBasesTrimmed -= s.length();
//        }
//        p.numRightBasesTrimmed = numBasesTrimmed;
        
        return p;
    }

    private synchronized PairedReadSegments nextFR() throws FileFormatException {
        PairedReadSegments p = nextFF();
        
        ArrayList<String> right = p.right;
        int numRightSegments = right.size();
        
        if (numRightSegments > 1) {
            Collections.reverse(right);
        }
        
        for (int i=0; i<numRightSegments; ++i) {
            right.set(i, reverseComplement(right.get(i)));
        }
        
        //p.right = right;
        
        return p;
    }

    private synchronized PairedReadSegments nextRF() throws FileFormatException {
        PairedReadSegments p = nextFF();

        ArrayList<String> left = p.left;
        int numLeftSegments = left.size();
        
        if (numLeftSegments > 1) {
            Collections.reverse(left);
        }
        
        for (int i=0; i<numLeftSegments; ++i) {
            left.set(i, reverseComplement(left.get(i)));
        }
        
        //p.left = left;
        
        return p;
    }

    private synchronized PairedReadSegments nextRR() throws FileFormatException {
        PairedReadSegments p = nextFF();
        
        ArrayList<String> left = p.left;
        int numLeftSegments = left.size();
        
        if (numLeftSegments > 1) {
            Collections.reverse(left);
        }
        
        for (int i=0; i<numLeftSegments; ++i) {
            left.set(i, reverseComplement(left.get(i)));
        }
        
        //p.left = left;
        
        ArrayList<String> right = p.right;
        int numRightSegments = right.size();
        
        if (numRightSegments > 1) {
            Collections.reverse(right);
        }
        
        for (int i=0; i<numRightSegments; ++i) {
            right.set(i, reverseComplement(right.get(i)));
        }
        
        //p.right = right;
        
        return p;
    }
    
    @Override
    public void close() throws IOException {
        leftReader.close();
        rightReader.close();
    }
}
