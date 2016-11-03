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
public class LargeBitBuffer extends AbstractLargeBitBuffer {
    private final long size;
    private final LargeByteBuffer backingByteBuffer;
    
    public LargeBitBuffer(long size) {
        this.size = size;
        long numBytes = size / Byte.SIZE;
        if (size % (long) Byte.SIZE > 0) {
            ++numBytes;
        }
        
        backingByteBuffer = new LargeByteBuffer(numBytes);
    }
    
    @Override
    public void set(long index) {
        long byteIndex = index / Byte.SIZE;
        byte b = backingByteBuffer.get(byteIndex);
        b |= (1 << (byte) (byteIndex % Byte.SIZE));
        backingByteBuffer.set(byteIndex, b);
    }

    @Override
    public boolean compareAndSwap(long index) {
        long byteIndex = index / Byte.SIZE;
        byte expected = backingByteBuffer.get(byteIndex);
        return backingByteBuffer.compareAndSwap(byteIndex, expected, (byte) (expected | (1 << (int) (index % Byte.SIZE))));    
    }
    
    @Override
    public boolean get(long index) {
        long byteIndex = index / Byte.SIZE;
        return (backingByteBuffer.get(byteIndex) & (1 << (byte) (byteIndex % Byte.SIZE))) != 0;
    }
    
    @Override
    public long size() {
        return size;
    }
    
    public void empty() {
        backingByteBuffer.empty();
    }
    
    @Override
    public long popCount() {
        return backingByteBuffer.bitPopCount();
    }
    
    @Override
    public void destroy() {
        backingByteBuffer.destroy();
    }
    
    @Override
    public void write(FileOutputStream out) throws IOException {
        backingByteBuffer.write(out);
    }
    
    @Override
    public void read(FileInputStream in) throws IOException {
        backingByteBuffer.read(in);
    }

    @Override
    public AbstractLargeByteBuffer getBackingByteBuffer() {
        return backingByteBuffer;
    }
}
