/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.buffer;

/**
 *
 * @author kmnip
 */
public abstract class AbstractLargeBitBuffer {
    public abstract void set(long index);
    public abstract boolean get(long index);
    public abstract long size();
    public abstract long popCount();
}
