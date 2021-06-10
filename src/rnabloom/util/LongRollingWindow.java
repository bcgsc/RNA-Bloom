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
    private int pos;
    private final int size;
    private int minIndex;
    
    public LongRollingWindow(long[] window) {
        this.window = window;
        this.size = window.length;
        this.index = size - 1; // assume the last value in the window is the last entry
        this.pos = index;
        updateMinIndex();
    }
    
    public void setIndex(int i, int pos) {
        this.index = i;
        this.pos = pos;
    }
    
    public void roll(long newVal) {
        ++pos;
        
        if (++index >= size) {
            index = 0;
        }
        
        window[index] = newVal;
        
        if (minIndex == index) {
            updateMinIndex();
        }
        else if (newVal < window[minIndex]) {
            minIndex = index;
        }
    }
    
    private void updateMinIndex() {
        minIndex = 0;
        long min = window[0];
        for (int i=1; i<size; ++i) {
            if (window[i] < min) {
                min = window[i];
                minIndex = i;
            }
        }
    }
    
    public long getMin() {        
        return window[minIndex];
    }
    
    public int getMinPos() {
        if (minIndex > index) {
            // window has wrapped around
            return pos - index - size + minIndex;
        }
        
        return pos - index + minIndex;
    }
}

