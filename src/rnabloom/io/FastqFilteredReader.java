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
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import static rnabloom.util.SeqUtils.getAveragePhred33Score;

/**
 *
 * @author Ka Ming Nip
 */
public class FastqFilteredReader extends FastqReader {
    private int minAvgBaseQual = 0;
        
    public FastqFilteredReader(String path, int minAvgBaseQual) throws IOException {
        super(path);
        this.minAvgBaseQual = minAvgBaseQual;
    }
    
    @Override
    public String next() throws FileFormatException {
        String line1, line3, seq, qual;
        
        try {
            synchronized(this) {
                line1 = itr.next();
                seq = itr.next();
                line3 = itr.next();
                qual = itr.next();
            }
        }
        catch (NoSuchElementException e) {
            return null;
        }
        catch (Exception e) {
            throw new FileFormatException("Error reading file");
        }
        
        if (line1.charAt(0) != '@') {
            throw new FileFormatException("Line 1 of FASTQ record is expected to start with '@'");
        }

        if (line3.charAt(0) != '+') {
            throw new FileFormatException("Line 3 of FASTQ record is expected to start with '+'");
        }
        
        if (minAvgBaseQual > getAveragePhred33Score(qual)) {
            seq = "";
        }
        
        return seq;
    }
    
    @Override
    public String[] nextWithName() throws FileFormatException {
        String line1, line3, name, seq, qual;
        
        try {
            synchronized(this) {
                line1 = itr.next();
                seq = itr.next();
                line3 = itr.next();
                qual = itr.next(); // line 4
            }
        }
        catch (NoSuchElementException e) {
            return null;
        }
        catch (Exception e) {
            throw new FileFormatException("Error reading file");
        }

        if (line3.charAt(0) != '+') {
            throw new FileFormatException("Line 3 of FASTQ record is expected to start with '+'");
        }
        
        Matcher m = RECORD_NAME_COMMENT_PATTERN.matcher(line1);
        if (m.matches()) {
            name = m.group(1);
            if (removeNameSuffix) {
                Matcher m2 = RECORD_NAME_PATTERN.matcher(name);
                if (m2.matches()) {
                    name = m2.group(1);
                }
            }
        }
        else {
            throw new FileFormatException("Line 1 of FASTQ record is expected to start with '@'");
        }
        
        if (minAvgBaseQual > getAveragePhred33Score(qual)) {
            seq = "";
        }
        
        return new String[]{name, seq, qual};
    }
    
    @Override
    public void nextWithoutName(FastqRecord fr) throws FileFormatException {        
        String line1, line3;
        
        try {
            synchronized(this) {
                line1 = itr.next();            
                fr.seq = itr.next();
                line3 = itr.next();
                fr.qual = itr.next();
            }
        }
        catch (NoSuchElementException e) {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            return;
        }
        catch (Exception e) {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            throw new FileFormatException("Error reading file");
        }

        if (line1.charAt(0) != '@') {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            throw new FileFormatException("Line 1 of FASTQ record is expected to start with '@'");
        }

        if (line3.charAt(0) != '+') {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            throw new FileFormatException("Line 3 of FASTQ record is expected to start with '+'");
        }
        
        if (minAvgBaseQual > getAveragePhred33Score(fr.qual)) {
            fr.seq = "";
            fr.seq = "";
        }
    }
    
    @Override
    public void nextWithName(FastqRecord fr) throws FileFormatException {
        String line1, line3;
        
        try {
            synchronized(this) {
                line1 = itr.next();
                fr.seq = itr.next();
                line3 = itr.next();            
                fr.qual = itr.next();
            }
        }
        catch (NoSuchElementException e) {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            return;
        }
        catch (Exception e) {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            throw new FileFormatException("Error reading file");
        }
        
        if (line3.charAt(0) != '+') {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            throw new FileFormatException("Line 3 of a FASTQ record is expected to start with '+'");
        }
        
        Matcher m = RECORD_NAME_COMMENT_PATTERN.matcher(line1);
        if (m.matches()) {
            fr.name = m.group(1);
            if (removeNameSuffix) {
                Matcher m2 = RECORD_NAME_PATTERN.matcher(fr.name);
                if (m2.matches()) {
                    fr.name = m2.group(1);
                }
            }
        }
        else {
            throw new FileFormatException("Line 1 of a FASTQ record is expected to start with '@'");
        }
        
        if (minAvgBaseQual > getAveragePhred33Score(fr.qual)) {
            fr.seq = "";
            fr.qual = "";
        }
    }
}
