/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

/**
 *
 * @author kmnip
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
