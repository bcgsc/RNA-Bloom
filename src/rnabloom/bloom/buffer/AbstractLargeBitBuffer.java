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
public abstract class AbstractLargeBitBuffer {
    public abstract void set(long index);
    public abstract void setCAS(long index);
    public abstract boolean get(long index);
    public abstract boolean getAndSet(long index);
    public abstract long size();
    public abstract long popCount();
    public abstract void empty();
    public abstract void destroy();
    public abstract void write(FileOutputStream out) throws IOException;
    public abstract void read(FileInputStream in) throws IOException;
    public abstract AbstractLargeByteBuffer getBackingByteBuffer();
}
