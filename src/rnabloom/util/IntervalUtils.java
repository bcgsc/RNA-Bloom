/* 
 * Copyright (C) 2021-present BC Cancer Genome Sciences Centre
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

import rnabloom.olc.Interval;

/**
 *
 * @author Ka Ming Nip
 */
public class IntervalUtils {
    
    public static boolean isContained(int start1, int end1, int start2, int end2) {
        return (start1 >= start2 && end1 <= end2) || (start2 >= start1 && end2 <= end1);
    }
    
    public static boolean isContained(Interval i1, Interval i2) {
        return isContained(i1.start, i1.end, i2.start, i2.end);
    }
    
    public static int getOverlap(int start1, int end1, int start2, int end2) {
        return Math.max(0, Math.min(end1, end2) - Math.max(start1, start2));
    }
    
    public static int getOverlap(Interval i1, Interval i2) {
        return getOverlap(i1.start, i1.end, i2.start, i2.end);
    }
        
    public static Interval merge(Interval i1, Interval i2) {
        if (getOverlap(i1, i2) > 0) {
            return new Interval(Math.min(i1.start, i2.start), Math.max(i1.end, i2.end));
        }
        return null;
    }
}
