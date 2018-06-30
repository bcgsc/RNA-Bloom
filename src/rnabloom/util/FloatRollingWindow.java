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

/**
 *
 * @author Ka Ming Nip
 */
public class FloatRollingWindow {
    private float[] window;
    private int index;
    private int size;
    
    public FloatRollingWindow(float[] window) {
        this.window = window;
        this.index = 0;
        this.size = window.length;
    }
    
    public void roll(float newVal) {
        if (++index >= size) {
            index = 0;
        }
        
        window[index] = newVal;
    }
    
    public float getMin() {
        float min = window[0];
        for (int i=1; i<size; ++i) {
            min = Math.min(min, window[i]);
        }
        
        return min;
    }
    
    public float getMax() {
        float max = window[0];
        for (int i=1; i<size; ++i) {
            max = Math.max(max, window[i]);
        }
        
        return max;
    }
}
