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
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import static rnabloom.io.Constants.BUFFER_SIZE;
import static rnabloom.util.FileUtils.getTextFileReader;

/**
 *
 * @author Ka Ming Nip
 */
public class PafReader {
    private final Iterator<String> itr;
    private final BufferedReader br;
    private long records = 0;
    
    public PafReader(String path) throws IOException {
        br = getTextFileReader(path);
        itr = br.lines().iterator();
    }
    
    public PafReader(InputStream in) {
        br = new BufferedReader(new InputStreamReader(in), BUFFER_SIZE);
        itr = br.lines().iterator();
    }
    
    public boolean hasNext() {
        return itr.hasNext();
    }
    
    public ExtendedPafRecord next() {
        ExtendedPafRecord r = new ExtendedPafRecord();
        r.update(itr.next());
        ++records;
        return r;
    }

    public void next(PafRecord record) {
        record.update(itr.next().trim().split("\t"));
        ++records;
    }
    
    public void next(ExtendedPafRecord record) {
        record.update(itr.next());
        ++records;
    }
    
    public void close() throws IOException {
        br.close();
    }
    
    public long getNumRecords() {
        return records;
    }
}
