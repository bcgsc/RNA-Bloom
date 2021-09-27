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
package rnabloom.util;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;

/**
 *
 * @author Ka Ming Nip
 */
public class Common {
    public static float convertToRoundedPercent(float f) {
        return roundToSigFigs(f * 100, 3);
    }
    
    public static float roundToSigFigs(float f, int sigFigs) {
        BigDecimal bd = new BigDecimal(f);
        bd = bd.round(new MathContext(sigFigs));
        return bd.floatValue();
    }
    
    public static float getMedian(float[] arr) {
        int len = arr.length;
        Arrays.sort(arr);
        int halfLen = len/2;
        if (len % 2 == 0) {
            return (arr[halfLen-1] + arr[halfLen])/2.0f;
        }
        
        return arr[halfLen];
    }
    
    public static float getMinium(float[] arr) {
        float min = Float.MAX_VALUE;
        for (float c : arr) {
            if (c < min) {
                min = c;
            }
        }
        return min;
    }
        
    public static float getMedian(Collection<Float> arr) {
        int len = arr.size();
        ArrayList<Float> a = new ArrayList<>(arr);
        Collections.sort(a);
        int halfLen = len/2;
        if (len % 2 == 0) {
            return (a.get(halfLen-1) + a.get(halfLen))/2.0f;
        }
        
        return a.get(halfLen);
    }

    public static float[] getMinMedMax(float[] a) {
        int len = a.length;
        Arrays.sort(a);
        int halfLen = len/2;
        if (len % 2 == 0) {
            return new float[]{a[0], (a[halfLen-1] + a[halfLen])/2.0f, a[len-1]};
        }
        
        return new float[]{a[0], a[halfLen], a[len-1]};
    }
    
    public static class Quartiles {
        public int min;
        public int q1;
        public int median;
        public int q3;
        public int max;
        
        @Override
        public String toString() {
            return "min:" + min + ", Q1:" + q1 + ", M:" + median + ", Q3:" + q3 + ", max:" + max;
        }
        
        public String toString(String sep) {
            return min + sep + q1 + sep + median + sep + q3 + sep + max;
        }
    }
    
    public static Quartiles getQuartiles(final ArrayList<Integer> arr) {
        int len = arr.size();
        
        Collections.sort(arr);
        
        Quartiles stats = new Quartiles();
        int halfLen = len/2;
        int q1Index = len/4;
        int q3Index = halfLen+q1Index;
                
        stats.min = arr.get(0);
        stats.max = arr.get(len-1);
        
        if (len % 2 == 0) {
            stats.median = (arr.get(halfLen-1) + arr.get(halfLen))/2;
        }
        else {
            stats.median = arr.get(halfLen);
        }
        
        if (len % 4 == 0) {
            stats.q1 = (arr.get(q1Index-1) + arr.get(q1Index))/2;
            stats.q3 = (arr.get(q3Index-1) + arr.get(q3Index))/2;
        }
        else {
            stats.q1 = arr.get(q1Index);
            stats.q3 = arr.get(q3Index);
        }
        
        return stats;
    }
    
    public static Quartiles getQuartiles(final int[] arr) {
        int len = arr.length;
        
        Arrays.sort(arr);
        
        Quartiles stats = new Quartiles();
        int halfLen = len/2;
        int q1Index = len/4;
        int q3Index = halfLen+q1Index;
                
        stats.min = arr[0];
        stats.max = arr[len-1];
        
        if (len % 2 == 0) {
            stats.median = (arr[halfLen-1] + arr[halfLen])/2;
        }
        else {
            stats.median = arr[halfLen];
        }
        
        if (len % 4 == 0) {
            stats.q1 = (arr[q1Index-1] + arr[q1Index])/2;
            stats.q3 = (arr[q3Index-1] + arr[q3Index])/2;
        }
        else {
            stats.q1 = arr[q1Index];
            stats.q3 = arr[q3Index];
        }
        
        return stats;
    }
}
