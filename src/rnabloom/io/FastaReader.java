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
import java.io.IOException;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static rnabloom.util.FileUtils.getTextFileReader;

/**
 *
 * @author Ka Ming Nip
 */
public class FastaReader implements FastxReaderInterface {
    private final static Pattern RECORD_NAME_PATTERN = Pattern.compile("([^\\s]+)/[12]");
    private final static Pattern RECORD_NAME_COMMENT_PATTERN = Pattern.compile("^>([^\\s]+)\\s*(.*)?$");
    private final Iterator<String> itr;
    private final BufferedReader br;
    private boolean removeNameSuffix = true;
    
    public FastaReader(String path) throws IOException {
        br = getTextFileReader(path);
        itr = br.lines().iterator();
    }

    public void setRemoveNameSuffix(boolean value) {
        this.removeNameSuffix = value;
    }
    
    public static boolean isCorrectFormat(String path) {
        try {
            // try to get the first FASTA record
            FastaReader reader = new FastaReader(path);
            String seq = reader.next();
            reader.close();
            if (seq == null) {
                return false;
            }
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

    private String header = null;
    
    @Override
    public synchronized String next() throws FileFormatException {
        if (!itr.hasNext()) {
            return null;
        }
        
        if (header == null) {
            header = itr.next().trim();
        }
        
        if (header.isEmpty()) {
            return null; //throw new FileFormatException("Empty FASTA header");
        }
        else if (header.charAt(0) != '>') {
            throw new FileFormatException("Incorrect FASTA header format");
        }
        
        StringBuilder builder = new StringBuilder();
        while (itr.hasNext()) {
            String line = itr.next().trim();
            if (line.isEmpty()) {
                header = null;
                break;
            }
            else if (line.charAt(0) == '>') {
                header = line;
                break;
            }
            else {
                builder.append(line);
            }
        }
        
        return builder.toString();
    }

    @Override
    public synchronized String[] nextWithName() throws FileFormatException {
        if (!itr.hasNext()) {
            return null;
        }
        
        if (header == null) {
            header = itr.next().trim();
        }
        
        if (header.isEmpty()) {
            return null; //throw new FileFormatException("Empty FASTA header");
        }
        
        String name = null;
        Matcher m = RECORD_NAME_COMMENT_PATTERN.matcher(header);
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
            throw new FileFormatException("Incorrect FASTA header format");
        }
        
        StringBuilder builder = new StringBuilder();
        while (itr.hasNext()) {
            String line = itr.next().trim();
            if (line.isEmpty()) {
                header = null;
                break;
            }
            else if (line.charAt(0) == '>') {
                header = line;
                break;
            }
            else {
                builder.append(line);
            }
        }
        
        return new String[]{name, builder.toString()};
    }
    
    public synchronized void nextWithName(FastaRecord fr) throws FileFormatException {
        if (!itr.hasNext()) {
            fr.name = null;
            fr.seq = null;
            return;
        }
        
        if (header == null) {
            header = itr.next().trim();
        }
        
        if (header.isEmpty()) {
            fr.name = null;
            fr.seq = null;
            return;
            //throw new FileFormatException("Empty FASTA header");
        }
        
        String name = null;
        Matcher m = RECORD_NAME_COMMENT_PATTERN.matcher(header);
        if (m.matches()) {
            name = m.group(1);
            Matcher m2 = RECORD_NAME_PATTERN.matcher(name);
            if (m2.matches()) {
                name = m2.group(1);
            }
        }
        else {
            throw new FileFormatException("Incorrect FASTA header format");
        }
        
        StringBuilder builder = new StringBuilder();
        while (itr.hasNext()) {
            String line = itr.next().trim();
            if (line.isEmpty()) {
                header = null;
                break;
            }
            else if (line.charAt(0) == '>') {
                header = line;
                break;
            }
            else {
                builder.append(line);
            }
        }
        
        fr.name = name;
        fr.seq = builder.toString();
    }
    
    public synchronized String[] nextWithComment() throws FileFormatException {
        if (!itr.hasNext()) {
            return null;
        }
        
        if (header == null) {
            header = itr.next().trim();
        }
        
        if (header.isEmpty()) {
            return null; //throw new FileFormatException("Empty FASTA header");
        }
    
        String name, comment = null;
        Matcher m = RECORD_NAME_COMMENT_PATTERN.matcher(header);
        if (m.matches()) {
            name = m.group(1);
            comment = m.group(2);
            if (removeNameSuffix) {
                Matcher m2 = RECORD_NAME_PATTERN.matcher(name);
                if (m2.matches()) {
                    name = m2.group(1);
                }
            }
        }
        else {
            throw new FileFormatException("Incorrect FASTA header format");
        }
        
        StringBuilder builder = new StringBuilder();
        while (itr.hasNext()) {
            String line = itr.next().trim();
            if (line.isEmpty()) {
                header = null;
                break;
            }
            else if (line.charAt(0) == '>') {
                header = line;
                break;
            }
            else {
                builder.append(line);
            }
        }
        
        return new String[]{name, comment, builder.toString()};
    }
    
    @Override
    public void close() throws IOException {
        br.close();
    }
    
    public static void main(String[] args) {
        //debug
    }
}
