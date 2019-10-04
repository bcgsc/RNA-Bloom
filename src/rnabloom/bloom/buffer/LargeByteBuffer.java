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
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.LongBuffer;

/**
 *
 * @author Ka Ming Nip
 */
public class LargeByteBuffer extends AbstractLargeByteBuffer {
    private final static int MAX_PARTITION_SIZE = Integer.MAX_VALUE;
    
    private final long size;
    private final ByteBuffer[] buffers;
    private int numPartitions;
    
    public LargeByteBuffer(long size) {
        this.size = size;
        this.numPartitions = (int) (size / MAX_PARTITION_SIZE);
        int remainder = (int) (size % MAX_PARTITION_SIZE);
        
        if (remainder > 0) {
            ++numPartitions;

            buffers = new ByteBuffer[numPartitions];
            
            int lastIndex = numPartitions-1;
            
            for (int i=0; i<lastIndex; ++i) {
                buffers[i] = ByteBuffer.allocateDirect(MAX_PARTITION_SIZE).order(ByteOrder.nativeOrder());
            }
            
            buffers[lastIndex] = ByteBuffer.allocateDirect(remainder).order(ByteOrder.nativeOrder());
        }
        else {
            buffers = new ByteBuffer[numPartitions];

            for (int i=0; i<numPartitions; ++i) {
                buffers[i] = ByteBuffer.allocateDirect(MAX_PARTITION_SIZE).order(ByteOrder.nativeOrder());
            }
        }
    }
    
    @Override
    public void set(long index, byte value) {
        buffers[(int) (index / MAX_PARTITION_SIZE)].put((int) (index % MAX_PARTITION_SIZE), value);
    }

    @Override
    public byte compareAndSwap(long index, byte expected, byte updated) {
        ByteBuffer bb = buffers[(int) (index / MAX_PARTITION_SIZE)];
        int bbIndex = (int) (index % MAX_PARTITION_SIZE);
        byte b = bb.get(bbIndex);
        if (expected == b) {
            bb.put(bbIndex, updated);
        }
        return b;
    }
    
    @Override
    public byte get(long index) {
        return buffers[(int) (index / MAX_PARTITION_SIZE)].get((int) (index % MAX_PARTITION_SIZE));
    }
    
    @Override
    public long size() {
        return size;
    }
    
    @Override
    public void empty() {
        for(ByteBuffer bb : buffers) {
            bb.clear();
        }
    }
    
    @Override
    public long popCount() {
        long count = 0;
        int cap;
        for(ByteBuffer bb : buffers) {
            cap = bb.capacity();
            for(int i=0; i<cap; ++i) {
                if (bb.get(i) > 0) {
                    ++count;
                }
            }
        }
        return count;
    }
    
    public long bitPopCount() {
        long count = 0;
        LongBuffer lb;
        int cap;
        for(ByteBuffer bb : buffers) {
            lb = bb.asLongBuffer();
            cap = lb.capacity();
            for(int i=0; i<cap; ++i) {
                count += Long.bitCount(lb.get(i));
            }
        }
        return count;
    }
    
    @Override
    public void destroy() {
        this.empty();
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
