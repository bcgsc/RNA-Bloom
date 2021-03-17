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

import rnabloom.util.SeqUtils;
import static rnabloom.util.SeqUtils.countBits;

/**
 *
 * @author Ka Ming Nip
 */
public class BufferComparator {
    public static boolean equivalentByteBuffers(AbstractLargeByteBuffer a, AbstractLargeByteBuffer b) {
        final long size = a.size();
        
        if (size != b.size()) {
            return false;
        }
        
        for (long i=0; i<size; ++i) {
            if (Byte.compare(a.get(i), b.get(i)) != 0) {
                return false;
            }
        }
        
        return true;
    }
    
    public static boolean equivalentBitBuffers(AbstractLargeBitBuffer a, AbstractLargeBitBuffer b) {        
        if (a.size() != b.size()) {
            return false;
        }
        
        return equivalentByteBuffers(a.getBackingByteBuffer(), b.getBackingByteBuffer());
    }
    
    public static long compareByteBuffers(AbstractLargeByteBuffer a, AbstractLargeByteBuffer b) {
        long size = a.size();
        
        if (size == b.size()) {
            long distance = 0;
            
            for (int i=0; i<size; ++i) {
                if (a.get(i) != b.get(i)) {
                    ++distance;
                }
            }
            
            return distance; // Hamming distance
        }
        
        return -1;
    }
    
    public static long compareBitBuffers(AbstractLargeBitBuffer a, AbstractLargeBitBuffer b) {
        if (a.size() == b.size()) {
            AbstractLargeByteBuffer aBytes = a.getBackingByteBuffer();
            AbstractLargeByteBuffer bBytes = b.getBackingByteBuffer();
            
            long size = aBytes.size();
            
            if (size == bBytes.size()) {
                long distance = 0;

                for (int i=0; i<size; ++i) {
                    distance += countBits((byte)(aBytes.get(i) ^ bBytes.get(i)));
                }

                return distance; // Hamming distance
            }
        }
        
        return -1;
    }
}
