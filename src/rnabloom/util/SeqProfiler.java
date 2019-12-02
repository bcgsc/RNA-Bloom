/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
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

import java.util.ArrayList;
import java.util.Collections;
import static rnabloom.util.SeqUtils.countN;
import static rnabloom.util.SeqUtils.countN;

/**
 *
 * @author Ka Ming Nip
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
