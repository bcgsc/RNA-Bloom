/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.olc;

import static rnabloom.util.IntervalUtils.getOverlap;

/**
 *
 * @author Ka Ming Nip
 */
public class ComparableInterval extends Interval implements Comparable<Interval> {
    public ComparableInterval(int start, int end) {
        super(start, end);
    }
    
    public boolean merge(Interval other) {
        return merge(other.start, other.end);
    }
    
    public boolean merge(int start2, int end2) {
        int overlap = getOverlap(start, end, start2, end2);
        if (overlap > 0) {
            start = Math.min(start, start2);
            end = Math.max(end, end2);
            
            return true;
        }
        return false;
    }

    @Override
    public int compareTo(Interval other) {
        int c = start - other.start;
        if (c == 0) {
            c = end - other.end;
        }
        return c;
    }
}
