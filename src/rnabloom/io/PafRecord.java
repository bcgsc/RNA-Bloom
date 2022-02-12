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

import rnabloom.olc.OverlapCoords;

/**
 *
 * @author Ka Ming Nip
 */
public class PafRecord extends OverlapCoords {
    public String qName = null, tName = null;
    public boolean reverseComplemented = false;
    public int qLen, tLen, numMatch, blockLen, qual = -1;
    
    public void update(String[] cols) {
        qName = cols[0];
        qLen = Integer.parseInt(cols[1]);
        qStart = Integer.parseInt(cols[2]);
        qEnd = Integer.parseInt(cols[3]);
        reverseComplemented = cols[4].equals("-");
        tName = cols[5];
        tLen = Integer.parseInt(cols[6]);
        tStart = Integer.parseInt(cols[7]);
        tEnd = Integer.parseInt(cols[8]);
        numMatch = Integer.parseInt(cols[9]);
        blockLen = Integer.parseInt(cols[10]);
        qual = Integer.parseInt(cols[11]);
    }
}
