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
package rnabloom.util;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author kmnip
 */
public class NTCardHistogram {
    public long numKmers = 0;
    public int covThreshold = 0;
    
    public NTCardHistogram(String path) throws FileNotFoundException, IOException {
        boolean f0Found = false;
        long lastCount = -1;
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.length() > 0) {
                String[] cols = line.split("\t");
                if (cols[0].equals("F0")) {
                    numKmers = Long.parseLong(cols[1]); 
                    f0Found = true;
                }
                else if (f0Found) {
                    long count = Long.parseLong(cols[1]);
                    if (lastCount <= 0) {
                        lastCount = count;
                    }
                    else if (3 * count > lastCount) {
                        covThreshold = Integer.parseInt(cols[0]) - 1;
                        break;
                    }
                    else {
                        lastCount = count;
                    }
                }
            }
        }
        br.close();
    }
    
}
