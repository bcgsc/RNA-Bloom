/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.buffer;

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
}
