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

import rnabloom.util.MiniFloat;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import rnabloom.bloom.buffer.UnsafeByteBuffer;
import rnabloom.bloom.buffer.AbstractLargeByteBuffer;
import rnabloom.bloom.buffer.LargeByteBuffer;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import rnabloom.bloom.buffer.BufferComparator;
import rnabloom.bloom.hash.HashFunction;

/**
 *
 * @author Ka Ming Nip
 */
public class CountingBloomFilter implements CountingBloomFilterInterface {
    protected AbstractLargeByteBuffer counts;
    protected int numHash;
    protected long size;
    protected HashFunction hashFunction;
    protected long popcount = -1;
    
    public CountingBloomFilter(long size, int numHash, HashFunction hashFunction) {
        this.size = size;
        try {
            //System.out.println("unsafe");
            this.counts = new UnsafeByteBuffer(size);
        }
        catch (NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.counts = new LargeByteBuffer(size);
        }
        this.numHash = numHash;
        this.hashFunction = hashFunction;
    }
    
    private static final String LABEL_SEPARATOR = ":";
    private static final String LABEL_SIZE = "size";
    private static final String LABEL_NUM_HASH = "numhash";
    private static final String LABEL_FPR = "fpr";
    
    public CountingBloomFilter(File desc, File bytes, HashFunction hashFunction) throws FileNotFoundException, IOException {        
        BufferedReader br = new BufferedReader(new FileReader(desc));
        String line;
        while ((line = br.readLine()) != null) {
            String[] entry = line.split(LABEL_SEPARATOR);
            String key = entry[0];
            String val = entry[1];
            switch(key) {
                case LABEL_SIZE:
                    size = Long.parseLong(val);
                    break;
                case LABEL_NUM_HASH:
                    numHash = Integer.parseInt(val);
                    break;
            }
        }
        br.close();
        
        this.hashFunction = hashFunction;
        
        try {
            //System.out.println("unsafe");
            this.counts = new UnsafeByteBuffer(size);
        }
        catch (NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.counts = new LargeByteBuffer(size);
        }
        
        FileInputStream fin = new FileInputStream(bytes);
        this.counts.read(fin);
        fin.close();
        
        /**@TODO Assert file size*/
    }
    
    private long getIndex(long hashVal) {
        // shift right to remove sign bit and modulus the size of buffer
        return (hashVal >>> 1) % size;
    }
    
    public void save(File desc, File bytes) throws IOException {
        FileWriter writer = new FileWriter(desc, false);
        
        writer.write(LABEL_SIZE + LABEL_SEPARATOR + this.size + "\n" +
                    LABEL_NUM_HASH + LABEL_SEPARATOR + this.numHash + "\n" +
                    LABEL_FPR + LABEL_SEPARATOR + this.getFPR() + "\n");
        writer.close();
        
        FileOutputStream out = new FileOutputStream(bytes, false);
        this.counts.write(out);
        out.close();
    }
        
    @Override
    public void increment(String key) {
        final long[] hashVals = new long[numHash];
        hashFunction.getHashValues(key, numHash, hashVals);
        increment(hashVals);
    }
        
    public void increment(final long[] hashVals) {
        // find the smallest count at all hash positions
        byte min = counts.get(getIndex(hashVals[0]));
        byte c;
        int h;
        for (h=1; h<numHash; ++h) {
            c = counts.get(getIndex(hashVals[h]));
            if (c < min) {
                min = c;
            }
            if (min == 0) {
                break;
            }
        }
        
        // increment the smallest count
        
        byte updated = MiniFloat.increment(min);
        
        if (updated != min) {
            // update min count only
            for (h=0; h<numHash; ++h) {
                counts.compareAndSwap(getIndex(hashVals[h]), min, updated);
            }
        }
    }

    @Override
    public float getCount(String key) {
        final long[] hashVals = new long[numHash];
        hashFunction.getHashValues(key, numHash, hashVals);
        return getCount(hashVals);
    }
    
    public float getCount(long hashVal) {
        return getCount(hashFunction.getHashValues(hashVal, numHash));
    }
    
    public float getCount(final long[] hashVals) {
        // find the smallest count
        byte min = counts.get(getIndex(hashVals[0]));
        byte c;
        for (int h=1; h<numHash; ++h) {
            c = counts.get(getIndex(hashVals[h]));
            if (c < min) {
                min = c;
            }
            if (min == 0) {
                return 0;
            }
        }
        
        // return the float value
        return MiniFloat.toFloat(min);
    }

    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
        popcount = counts.popCount();
        return (float) pow((double)(popcount) / (double)(size), numHash);
    }

    public static long getExpectedSize(long expNumElements, float fpr, int numHash) {
        double r = (double) (-numHash) / log(1 - exp(log(fpr) / (double) numHash));
        return (long) Math.ceil(expNumElements * r);
    }
    
    public long getOptimalSize(float fpr) {
        if (popcount > 0) {
            double r = (double) (-numHash) / log(1 - exp(log(fpr) / (double) numHash));
            return (long) Math.ceil(popcount * r);
        }
        else {
            return size;
        }
    }
    
    public long getPopCount() {
        return popcount;
    }
    
//    public long updatePopcount() {
//        popcount = counts.popCount();
//        return popcount;
//    }
    
//    public long getSizeNeeded(float fpr) {
//        if (popcount < 0) {
//            popcount = updatePopcount();
//        }
//        
//        return (long) Math.ceil(- popcount *log(fpr) / pow(log(2), 2));
//    }
    
    public int getNumHash() {
        return numHash;
    }
    
//    public float getOptimalNumHash() {
//        if (popcount < 0) {
//            popcount = updatePopcount();
//        }
//        
//        return (float) (size / (float) popcount * log(2));
//    }
    
    public void empty() {
        if (this.counts != null) {
            this.counts.empty();
        }
    }
    
    public void destroy() {
        if (this.counts != null) {
            this.counts.destroy();
            this.counts = null;
        }
    }
    
    public boolean equivalent(CountingBloomFilter bf) {
        return this.size == bf.size &&
                this.numHash == bf.numHash &&
                BufferComparator.equivalentByteBuffers(counts, bf.counts);
    }
}
