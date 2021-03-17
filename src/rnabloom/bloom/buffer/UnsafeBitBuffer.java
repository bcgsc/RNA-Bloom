/* 
 * Copyright (C) 2018-present BC Cancer Genome Sciences Centre
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

/**
 *
 * @author Ka Ming Nip
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
        
    private static byte getBitMask(long i) {
        return (byte) (1 << (int) (i % Byte.SIZE));
    }
    
    private static long getByteBufferIndex(long i) {
        return i / Byte.SIZE;
    }
    
    @Override
    public void set(long index) {
        backingByteBuffer.or(getByteBufferIndex(index), getBitMask(index));
    }
    
    @Override
    public void setCAS(long index) {
        long byteIndex = getByteBufferIndex(index);
        byte expected = backingByteBuffer.get(byteIndex);
        byte mask = getBitMask(index);
        while (true) {
            byte b = backingByteBuffer.compareAndSwap(byteIndex, expected, (byte) (expected | mask));
            if (b == expected) {
                return;
            }
            expected = b;
        }
    }

    @Override
    public boolean get(long index) {
        return (backingByteBuffer.get(getByteBufferIndex(index)) & getBitMask(index)) != 0;
    }

    @Override
    public boolean getAndSet(long index) {
        return backingByteBuffer.compareAndOr(getByteBufferIndex(index), getBitMask(index));
    }
    
    @Override
    public long size() {
        return size;
    }
    
    @Override
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
