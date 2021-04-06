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

/**
 *
 * @author Ka Ming Nip
 */
public class LongRollingWindow {
    private final long[] window;
    private int index;
    private final int size;
    
    public LongRollingWindow(long[] window) {
        this.window = window;
        this.index = 0;
        this.size = window.length;
    }
    
    public void roll(long newVal) {
        if (++index >= size) {
            index = 0;
        }
        
        window[index] = newVal;
    }
    
    public long getMin() {
        long min = window[0];
        for (int i=1; i<size; ++i) {
            min = Math.min(min, window[i]);
        }
        
        return min;
    }
    
    public long getMax() {
        long max = window[0];
        for (int i=1; i<size; ++i) {
            max = Math.max(max, window[i]);
        }
        
        return max;
    }
}

