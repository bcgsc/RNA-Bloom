/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.buffer;

import java.lang.reflect.Field;
import sun.misc.Unsafe;

/**
 *
 * @author kmnip
 */
public class UnsafeByteBuffer extends AbstractLargeByteBuffer {
    private final long start;
    private final long size;
    private final Unsafe unsafe;
    
    private static Unsafe getMyUnsafe() throws NoSuchFieldException, IllegalArgumentException, IllegalAccessException {
        Field theUnsafe = Unsafe.class.getDeclaredField("theUnsafe");
        theUnsafe.setAccessible(true);
        return (Unsafe) theUnsafe.get(null);
    }
    
    public UnsafeByteBuffer(long size) throws NoSuchFieldException, IllegalArgumentException, IllegalAccessException {
        unsafe = getMyUnsafe();
        this.start = unsafe.allocateMemory(size);
        this.size = size;
        this.empty();
    }

    @Override
    public void set(long index, byte value) {
        unsafe.putByte(start + index, value);
    }

    @Override
    public byte get(long index) {
        return unsafe.getByte(start + index);
    }
    
    @Override
    public long size() {
        return size;
    }
    
    public void empty() {
        unsafe.setMemory(start, size, (byte) 0);
    }
    
    @Override
    public long popCount() {
        long count = 0;
        for (long i=0; i<size; ++i) {
            if (get(i) > 0) {
                ++count;
            }
        }
        return count;
    }
    
    public long bitPopCount() {
        long count = 0;
        for (long i=0; i<size; ++i) {
            byte b = get(i);
            while (b != 0) {
                count += (b & 1);
                b = (byte) (b >> 1);
            }
        }
        return count;
    }
    
    public void destroy() {
        unsafe.freeMemory(start);
    }
}
