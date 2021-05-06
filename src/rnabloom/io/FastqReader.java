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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.NoSuchElementException;
//import java.util.function.Supplier;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import static rnabloom.io.Constants.BUFFER_SIZE;
import static rnabloom.io.Constants.GZIP_EXTENSION;

/**
 *
 * @author Ka Ming Nip
 */
public final class FastqReader implements FastxReaderInterface {
    private final static Pattern RECORD_NAME_PATTERN = Pattern.compile("([^\\s]+)/[12]");
    private final static Pattern RECORD_NAME_COMMENT_PATTERN = Pattern.compile("^@([^\\s]+)\\s*(.*)?$");
    private final BufferedReader br;
    private final Iterator<String> itr;
    
    public FastqReader(String path) throws IOException {        
        if (path.toLowerCase().endsWith(GZIP_EXTENSION)) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))), BUFFER_SIZE);
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(path)), BUFFER_SIZE);
        }
        itr = br.lines().iterator();
    }

    public static boolean isCorrectFormat(String path) {
        try {
            // try to get the first FASTQ record
            FastqReader reader = new FastqReader(path);
            reader.nextWithoutName(new FastqRecord());
            reader.close();
        }
        catch (Exception e) {
            return false;
        }
        
        return true;
    }
    
    @Override
    public synchronized boolean hasNext() {
        return itr.hasNext();
    }
    
    @Override
    public synchronized String next() throws FileFormatException {
        if (!itr.hasNext()) {
            return null;
        }
        
        String line1, line3, seq;
        
        synchronized(this) {
            line1 = itr.next();
            seq = itr.next();
            line3 = itr.next();
            itr.next();
        }
        
        if (line1.charAt(0) != '@') {
            throw new FileFormatException("Line 1 of FASTQ record is expected to start with '@'");
        }

        if (line3.charAt(0) != '+') {
            throw new FileFormatException("Line 3 of FASTQ record is expected to start with '+'");
        }
        
        return seq;
    }
    
    @Override
    public synchronized String[] nextWithName() throws FileFormatException {
        if (!itr.hasNext()) {
            return null;
        }
        
        String line1, line3, name, seq;
        
        try {
            synchronized(this) {
                line1 = itr.next();
                seq = itr.next();
                line3 = itr.next();
                itr.next(); // line 4
            }
        }
        catch (NoSuchElementException e) {
            return null;
        }

        if (line3.charAt(0) != '+') {
            throw new FileFormatException("Line 3 of FASTQ record is expected to start with '+'");
        }
        
        Matcher m = RECORD_NAME_COMMENT_PATTERN.matcher(line1);
        if (m.matches()) {
            name = m.group(1);
            Matcher m2 = RECORD_NAME_PATTERN.matcher(name);
            if (m2.matches()) {
                name = m2.group(1);
            }
        }
        else {
            throw new FileFormatException("Line 1 of FASTQ record is expected to start with '@'");
        }
        
        return new String[]{name, seq};
    }
    
    public synchronized void nextWithoutName(FastqRecord fr) throws FileFormatException {
        if (!itr.hasNext()) {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            return;
        }
        
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

        if (line1.charAt(0) != '@') {
            throw new FileFormatException("Line 1 of FASTQ record is expected to start with '@'");
        }

        if (line3.charAt(0) != '+') {
            throw new FileFormatException("Line 3 of FASTQ record is expected to start with '+'");
        }
        
    }
    
    public synchronized void nextWithName(FastqRecord fr) throws FileFormatException {
        if (!itr.hasNext()) {
            fr.name = null;
            fr.qual = null;
            fr.seq = null;
            return;
        }
        
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
        
        if (line3.charAt(0) != '+') {
            throw new FileFormatException("Line 3 of a FASTQ record is expected to start with '+'");
        }
        
        Matcher m = RECORD_NAME_COMMENT_PATTERN.matcher(line1);
        if (m.matches()) {
            fr.name = m.group(1);
            Matcher m2 = RECORD_NAME_PATTERN.matcher(fr.name);
            if (m2.matches()) {
                fr.name = m2.group(1);
            }
        }
        else {
            throw new FileFormatException("Line 1 of a FASTQ record is expected to start with '@'");
        }
    }
        
    @Override
    public void close() throws IOException {
        br.close();
    }
    
    public static void main(String[] args) {
        //debug
    }
}
