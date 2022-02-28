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

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import rnabloom.olc.ComparableInterval;
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
    
    public static boolean isForwardDoveTail(int start1, int end1, int start2, int end2) {
        return start1 < start2 && end1 < end2 && end1 > start2;
    }
    
    public static boolean isDoveTail(int start1, int end1, int start2, int end2) {
        return (start1 < start2 && end1 < end2 && end1 > start2) || 
                (start2 < start1 && end2 < end1 && end2 > start1);
    }
    
    public static boolean isDoveTail(Interval i1, Interval i2) {
        return isDoveTail(i1.start, i1.end, i2.start, i2.end);
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
    
    public static ArrayDeque<ComparableInterval> mergeIntervals(Collection<ComparableInterval> intervals) {
        ArrayList<ComparableInterval> intervalsList = new ArrayList<>(intervals);
        Collections.sort(intervalsList);
        
        ArrayDeque<ComparableInterval> newList = new ArrayDeque<>();
        newList.add(intervalsList.get(0));
        
        int numIntervals = intervalsList.size();
        for (int i=1; i<numIntervals; ++i) {
            ComparableInterval interval = intervalsList.get(i);
            if (!newList.getLast().merge(interval)) {
                newList.add(interval);
            }
        }
        
        return newList;
    }
}
