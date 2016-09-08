/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.buffer;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Field;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Scanner;
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
    
    @Override
    public void destroy() {
        unsafe.freeMemory(start);
    }
    
    private final static int TMP_BUFF_SIZE = 50000000;
    
    @Override
    public void write(FileOutputStream out) throws IOException {        
        byte[] buffer = new byte[TMP_BUFF_SIZE];
        long i = 0;
        
        while (i < size) {
            if (i+TMP_BUFF_SIZE > size) {
                int bufSize = (int) (size-i);
                buffer = new byte[bufSize];
                
                for (int bi=0; bi<bufSize; ++bi) {
                    buffer[bi] = this.get(i++);
                }                
            }
            else {
                for (int bi=0; bi<TMP_BUFF_SIZE; ++bi) {
                    buffer[bi] = this.get(i++);
                }
            }
            
            out.write(buffer);
        }
    }
    
    @Override
    public void read(FileInputStream in) throws IOException {
        byte[] buffer = new byte[TMP_BUFF_SIZE];
        long i = 0;
        
        while (i < size) {
            if (i+TMP_BUFF_SIZE > size) {
                buffer = new byte[(int) (size-i)];
            }
            
            in.read(buffer);
            
            for (byte b : buffer) {
                this.set(i++, b);
            }
        }
    }
}
