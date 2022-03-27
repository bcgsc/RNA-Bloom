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

/**
 *
 * @author Ka Ming Nip
 */
public class ExtendedPafRecord extends PafRecord {
    public String cigar = null; // CIGAR string, eg. "129M" in "cg:Z:129M"
    public int nm = -1; // Total number of mismatches and gaps in the alignment, eg. 2 in "NM:i:2"
    public boolean isPrimary = false;
    private String line = null;
    
    public void update(String line) {
        this.line = line;
        String[] cols = line.trim().split("\t");
        super.update(cols);
        
        // reset to default values then update
        nm = -1;
        cigar = null;
        isPrimary = false;
        
        int numCols = cols.length;
        for (int i=12; i<numCols; ++i) {
            String[] triplet = cols[i].split(":");
            String tag = triplet[0];
            switch(tag) {
                case "NM":
                    nm = Integer.parseInt(triplet[2]);
                    break;
                case "cg":
                    cigar = triplet[2];
                    break;
                case "tp":
                    isPrimary = triplet[2].equals("P");
                    break;
            }
        }
    }
    
    @Override
    public String toString() {
        return line;
    }
}
