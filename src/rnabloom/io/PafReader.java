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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import static rnabloom.io.Constants.BUFFER_SIZE;
import static rnabloom.io.Constants.GZIP_EXTENSION;

/**
 *
 * @author kmnip
 */
public class PafReader {
    private final Iterator<String> itr;
    private final BufferedReader br;
    
    public PafReader(String path) throws IOException {
        if (path.toLowerCase().endsWith(GZIP_EXTENSION)) {
            br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(path))), BUFFER_SIZE);
        }
        else {
            br = new BufferedReader(new InputStreamReader(new FileInputStream(path)), BUFFER_SIZE);
        }
        itr = br.lines().iterator();
    }
    
    public PafReader(InputStream in) {
        br = new BufferedReader(new InputStreamReader(in));
        itr = br.lines().iterator();
    }
    
    public boolean hasNext() {
        return itr.hasNext();
    }
    
    public ExtendedPafRecord next() {
        ExtendedPafRecord r = new ExtendedPafRecord();
        r.update(itr.next().trim().split("\t"));
        return r;
    }
    
    public void next(ExtendedPafRecord record) {
        String[] cols = itr.next().trim().split("\t");
        record.update(cols);
    }
    
    public void close() throws IOException {
        br.close();
    }
}
