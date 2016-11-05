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
public class UnsafeBitBuffer extends AbstractLargeBitBuffer {
    
    private final long size;
    private final UnsafeByteBuffer backingByteBuffer;
    
    public UnsafeBitBuffer(long size) throws NoSuchFieldException, IllegalArgumentException, IllegalAccessException {
        this.size = size;
        long numBytes = size / Byte.SIZE;
        if (size % Byte.SIZE > 0) {
            ++numBytes;
        }
        
        backingByteBuffer = new UnsafeByteBuffer(numBytes);
    }
    
    @Override
    public void set(long index) {
        backingByteBuffer.or(index / Byte.SIZE, (byte) (1 << (int) (index % Byte.SIZE)));
    }
    
    @Override
    public void setCAS(long index) {
        long byteIndex = index / Byte.SIZE;
        byte expected = backingByteBuffer.get(byteIndex);
        while (!backingByteBuffer.compareAndSwap(byteIndex, expected, (byte) (expected | (1 << (int) (index % Byte.SIZE))))) {
            expected = backingByteBuffer.get(byteIndex);
        }
    }

    @Override
    public boolean get(long index) {
        long byteIndex = index / Byte.SIZE;
        return (backingByteBuffer.get(byteIndex) & (1 << (int) (index % Byte.SIZE))) != 0;
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
