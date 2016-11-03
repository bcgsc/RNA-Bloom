/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.buffer;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 *
 * @author kmnip
 */
public abstract class AbstractLargeByteBuffer {
    public abstract void set(long index, byte value);
    public abstract boolean compareAndSwap(long index, byte expected, byte updated);
    public abstract byte get(long index);
    public abstract long size();
    public abstract long popCount();
    public abstract void destroy();
    public abstract void write(FileOutputStream out) throws IOException;
    public abstract void read(FileInputStream in) throws IOException;
}
