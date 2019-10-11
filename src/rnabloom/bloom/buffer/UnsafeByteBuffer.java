/* 
 * Copyright (C) 2018 BC Cancer Genome Sciences Centre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package rnabloom.bloom.buffer;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.lang.reflect.Field;
import sun.misc.Unsafe;

/**
 *
 * @author Ka Ming Nip
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

    public void or(long index, byte value) {
        long i = start + index;
        unsafe.putByte(i, (byte) (unsafe.getByte(i) | value));
    }
    
    public boolean compareAndOr(long index, byte value) {
        long i = start + index;
        byte b = unsafe.getByte(i);
        byte newB = (byte) (b | value);
        if (b != newB) {
            unsafe.putByte(i, newB);
            return false;
        }
        
        return true;
    }
    
    /*
    // rotate "v" to the left by "s" positions
    public static int rol(final int v, final int s) {
        return (v << s) | (v >>> (32 - s));
    }

    // rotate "v" to the right by "s" positions
    public static int ror(final int v, final int s) {
        return (v >>> s) | (v << (32 - s));
    }
    
    // rotate "v" to the left by "s" positions
    public static long rol(final long v, final int s) {
        return (v << s) | (v >>> (64 - s));
    }

    // rotate "v" to the right by "s" positions
    public static long ror(final long v, final int s) {
        return (v >>> s) | (v << (64 - s));
    }
    */

    @Override
    public byte compareAndSwap(long index, byte expected, byte updated) {
        long i = start + index;
        
        byte b = unsafe.getByte(i);
        if (expected == b) {
            unsafe.putByte(i, updated);
        }
                
        return b;
    }

    @Override
    public byte get(long index) {
        return unsafe.getByte(start + index);
    }
    
    @Override
    public long size() {
        return size;
    }
    
    @Override
    public void empty() {
        unsafe.setMemory(start, size, (byte) 0);
    }
    
    @Override
    public long popCount() {
        long count = 0;
        for (long i=0; i<size; ++i) {
            if (unsafe.getByte(start + i) != 0) {
                ++count;
            }
        }
        return count;
    }
        
    public long bitPopCount() {
        long count = 0;
        
        long numLongs = (long) Math.floor(size / 8f);
        for (long i=0; i<numLongs; ++i) {
            count += Long.bitCount(unsafe.getLong(start + (i * 8L)));
        }
        
        for (long i=numLongs * 8L; i<size; ++i) {
            byte b = get(i);

            b = (byte) ((byte) (b & 0b01010101) + (byte)(((byte) (b >>> 1)) & 0b01010101));
            b = (byte) ((byte) (b & 0b00110011) + (byte)(((byte) (b >>> 2)) & 0b00110011));
            b = (byte) ((byte) (b & 0b00001111) + (byte)(((byte) (b >>> 4)) & 0b00001111));
            
            count += b;
        }
        
        return count;
    }
    
    @Override
    public void destroy() {
        unsafe.freeMemory(start);
    }
    
    private final static int TMP_BUFF_SIZE = 1000000000; // 1 GB
    
    @Override
    public void write(FileOutputStream out) throws IOException {        
        byte[] buffer = size>TMP_BUFF_SIZE ? new byte[TMP_BUFF_SIZE] : new byte[(int)size];
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
        
        out.flush();
    }
    
    @Override
    public void read(FileInputStream in) throws IOException {
        byte[] buffer = size>TMP_BUFF_SIZE ? new byte[TMP_BUFF_SIZE] : new byte[(int)size];
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
