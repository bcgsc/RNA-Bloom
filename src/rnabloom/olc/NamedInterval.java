/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.olc;

import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author Ka Ming Nip
 */
public class NamedInterval extends Interval implements Comparable<Interval> {
    public String name;
    
    public NamedInterval(String name, int start, int end) {
        super(start, end);
        this.name = name;
    }

    @Override
    public int compareTo(Interval o) {
        return (end - start) - (o.end - o.start);
    }
    
    public static void main(String[] args) {
        NamedInterval i1 = new NamedInterval("1", 0, 10);
        NamedInterval i2 = new NamedInterval("2", 1, 9);
        NamedInterval i3 = new NamedInterval("3", 2, 8);
        ArrayList<NamedInterval> list = new ArrayList<>();
        list.add(i3);
        list.add(i2);
        list.add(i1);
        Collections.sort(list, Collections.reverseOrder());
        System.out.println(list);
    }
}
