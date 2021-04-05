/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.olc;

/**
 *
 * @author Ka Ming Nip
 */
public class ComparableInterval extends Interval implements Comparable<Interval> {
    public ComparableInterval(int start, int end) {
        super(start, end);
    }

    public static Interval merge(Interval a, Interval b) {
        if (a.start >= b.start && a.start <= b.end) {
            return new Interval(b.start, Math.max(a.end, b.end));
        }
        else if (a.end >= b.start && a.end <= b.end) {
            return new Interval(Math.min(a.start, b.start), b.end);
        }
        else if (b.start >= a.start && b.start <= a.end) {
            return new Interval(a.start, Math.max(a.end, b.end));
        }
        else if (b.end >= a.start && b.end <= a.end) {
            return new Interval(Math.min(a.start, b.start), b.end);
        }
        return null; // not merged
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
