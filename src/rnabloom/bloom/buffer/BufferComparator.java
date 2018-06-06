/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.buffer;

import rnabloom.util.SeqUtils;
import static rnabloom.util.SeqUtils.countBits;

/**
 *
 * @author kmnip
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
