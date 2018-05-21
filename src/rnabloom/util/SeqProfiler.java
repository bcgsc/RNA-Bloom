/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.util.ArrayList;
import java.util.Collections;
import static rnabloom.util.SeqUtils.countN;

/**
 *
 * @author kmnip
 */
public class SeqProfiler {
    private int size = 0;
    
    private ArrayList<Integer> lengths = new ArrayList<>();
    private ArrayList<Integer> a = new ArrayList<>();
    private ArrayList<Integer> c = new ArrayList<>();
    private ArrayList<Integer> g = new ArrayList<>();
    private ArrayList<Integer> t = new ArrayList<>();
//    private ArrayList<Integer> gc = new ArrayList<>();
    
    public SeqProfiler() {
    }
    
    private final int[] tempResult = new int[4];
    
    public void add(String seq) {
        lengths.add(seq.length());
        
        countN(seq, tempResult);
        
        a.add(tempResult[0]);
        c.add(tempResult[1]);
        g.add(tempResult[2]);
        t.add(tempResult[3]);
        
        ++size;
    }
    
    public int getSize() {
        return size;
    }

    public int[] evalLen() {
        return getStats(lengths);
    }
    
    public int[] evalA() {
        return getStats(a);
    }

    public int[] evalC() {
        return getStats(c);
    }

    public int[] evalG() {
        return getStats(g);
    }
    
    public int[] evalT() {
        return getStats(t);
    }
    
//    public int[] evalGC() {
//        return getStats(gc);
//    }
    
    private int[] getStats(ArrayList<Integer> list) {
        if (size <= 0) {
            return null;
        }
        
        int halfLen = size/2;
        int q1Index = size/4;
        int q3Index = halfLen+q1Index;
        
        int[] result = new int[5];
        Collections.sort(list);
        
        result[0] = list.get(0);
        result[4] = list.get(size-1);
        
        if (size % 2 == 0) {
            result[2] = (list.get(halfLen-1) + list.get(halfLen))/2;
        }
        else {
            result[2] = list.get(halfLen);
        }
        
        if (size % 4 == 0) {
            result[1] = (list.get(q1Index-1) + list.get(q1Index))/2;
            result[3] = (list.get(q3Index-1) + list.get(q3Index))/2;
        }
        else {
            result[1] = list.get(q1Index);
            result[3] = list.get(q3Index);
        }
        
        return result;
    }
}
