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
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.LongBuffer;

/**
 *
 * @author kmnip
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
    public byte get(long index) {
        return buffers[(int) (index / MAX_PARTITION_SIZE)].get((int) (index % MAX_PARTITION_SIZE));
    }
    
    @Override
    public long size() {
        return size;
    }
    
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
        
        out.flush();
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
