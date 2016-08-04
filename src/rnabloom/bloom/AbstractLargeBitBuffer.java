/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

/**
 *
 * @author kmnip
 */
abstract class AbstractLargeBitBuffer {
    abstract void set(long index);
    abstract boolean get(long index);
    abstract long size();
    abstract long popCount();
}
