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
package rnabloom.util;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import static rnabloom.util.FileUtils.getTextFileReader;

/**
 *
 * @author Ka Ming Nip
 */
public class NTCardHistogram {
    public static final int MAX_MULTIPLICITY = 65535;
    public long numKmers = 0;
    public long numUniqueKmers = 0;
    public long numUniqueOverrepresentedKmers = 0;
    public long[] counts = new long[MAX_MULTIPLICITY];
    
    public NTCardHistogram(String path) throws FileNotFoundException, IOException {
        boolean f0Found = false;
        BufferedReader br = getTextFileReader(path);
        String line;
        long sum = 0;
        while ((line = br.readLine()) != null) {
            if (line.length() > 0) {
                String[] cols = line.split("\t");
                if (f0Found) {
                    long count = Long.parseLong(cols[1]);
                    counts[Integer.parseInt(cols[0]) - 1] = count;
                    sum += count;
                }
                else if (cols[0].equals("F1")) {
                    numKmers = Long.parseLong(cols[1]);
                }
                else if (cols[0].equals("F0")) {
                    numUniqueKmers = Long.parseLong(cols[1]); 
                    f0Found = true;
                    if ((numUniqueKmers < 1 && numKmers > 0) ||
                            numKmers < numUniqueKmers) {
                        break;
                    }
                }
            }
        }
        
        numUniqueOverrepresentedKmers = numUniqueKmers - sum;
        br.close();
    }
    
    public long getNumSingletons() {
        return counts[0];
    }
    
    public int getMinCovThreshold(int multiplier) {
        for (int i=1; i<MAX_MULTIPLICITY; ++i) {
            if (multiplier * counts[i] > counts[i-1]) {
                return i;
            }
        }
        
        return 0;
    }
    
    public int getMaxCovThreshold(double fraction) {
        long numKmers = Math.round(fraction * numUniqueKmers);
        long sum = numUniqueOverrepresentedKmers;
        if (sum >= numKmers) {
            return MAX_MULTIPLICITY+1;
        }
        
        for (int i=MAX_MULTIPLICITY-1; i>=0; --i) {
            sum += counts[i];
            if (sum >= numKmers) {
                return i+1;
            }
        }
        
        return MAX_MULTIPLICITY + 1;
    }
}
