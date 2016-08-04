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
public class LargeBitBuffer {
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
    
    public void set(long index) {
        long byteIndex = index / Byte.SIZE;
        byte b = backingByteBuffer.get(byteIndex);
        b |= (1 << (byte) (byteIndex % Byte.SIZE));
        backingByteBuffer.set(byteIndex, b);
    }

    public boolean get(long index) {
        long byteIndex = index / Byte.SIZE;
        return (backingByteBuffer.get(byteIndex) & (1 << (byte) (byteIndex % Byte.SIZE))) != 0;
    }
    
    public long size() {
        return size;
    }
    
    public void empty() {
        backingByteBuffer.empty();
    }
    
    public long popCount() {
        return backingByteBuffer.bitPopCount();
    }
}
