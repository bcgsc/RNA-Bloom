/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.olc;

/**
 *
 * @author gengar
 */
public class ComparableInterval extends Interval implements Comparable<Interval> {
    public ComparableInterval(int start, int end) {
        super(start, end);
    }

    public boolean merge(Interval other) {
        if (start >= other.start && start <= other.end) {
            start = other.start;
            end = Math.max(end, other.end);
            return true;
        }
        else if (end >= other.start && end <= other.end) {
            start = Math.min(start, other.start);
            end = other.end;
            return true;
        }
        else if (other.start >= start && other.start <= end) {
            end = Math.max(end, other.end);
            return true;
        }
        else if (other.end >= start && other.end <= end) {
            start = Math.min(start, other.start);
            return true;
        }
        return false; // not merged
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
