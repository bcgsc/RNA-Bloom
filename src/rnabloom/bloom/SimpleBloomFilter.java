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
package rnabloom.bloom;

import rnabloom.bloom.buffer.AbstractLargeBitBuffer;
import rnabloom.bloom.buffer.LargeBitBuffer;
import rnabloom.bloom.buffer.UnsafeBitBuffer;

/**
 *
 * @author Ka Ming Nip
 */
public class SimpleBloomFilter implements BloomFilterInterface {
    protected AbstractLargeBitBuffer bitArray;
    protected long size;
    protected long popcount = -1;
    
    public SimpleBloomFilter(long size) {
        this.size = size;
        try {
            this.bitArray = new UnsafeBitBuffer(size);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArray = new LargeBitBuffer(size);
        }
    }
    
    protected long getIndex(String key) {
        long hashCode = (long) key.hashCode();
        return (hashCode + (long) Integer.MAX_VALUE + 1L) % size;
    }
    
    @Override
    public void add(String key) {
        bitArray.set(getIndex(key));
    }

    @Override
    public boolean lookup(String key) {
        return bitArray.get(getIndex(key));
    }
    
    public boolean lookupAndAdd(String key) {
        return bitArray.getAndSet(getIndex(key));
    }

    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
        popcount = bitArray.popCount();
        return (float) ((double) popcount / (double) size);
    }
    
    public void empty() {
        if (this.bitArray != null) {
            this.bitArray.empty();
        }
    }
    
    public void destroy() {
        if (this.bitArray != null) {
            this.bitArray.destroy();
            this.bitArray = null;
        }
    }
}
