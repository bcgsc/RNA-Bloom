/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

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
//import static java.lang.Math.exp;
//import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.random;
import static java.lang.Math.scalb;
import rnabloom.bloom.buffer.BufferComparator;
import rnabloom.bloom.hash.HashFunction2;

/**
 *
 * @author kmnip
 */
public class CountingBloomFilter implements CountingBloomFilterInterface {
    protected AbstractLargeByteBuffer counts;
    protected int numHash;
    protected long size;
    protected HashFunction2 hashFunction;
    protected long popcount = -1;
        
    private static final byte MANTISSA = 3;
    private static final byte MANTI_MASK = 0xFF >> (8 - MANTISSA);
    private static final byte ADD_MASK = 0x80 >> (7 - MANTISSA);
    
    public CountingBloomFilter(long size, int numHash, HashFunction2 hashFunction) {
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
    
    public CountingBloomFilter(File desc, File bytes, HashFunction2 hashFunction) throws FileNotFoundException, IOException {        
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
        if (min <= MANTI_MASK ||
                (min < Byte.MAX_VALUE &&
                (int) (random() * Integer.MAX_VALUE) % (1 << ((min >> MANTISSA) - 1)) == 0)) {
            byte updated = (byte) (min + 1);
            
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
        if (min <= MANTI_MASK) {
            return (float) min;
        }
        
        return scalb((min & MANTI_MASK) | ADD_MASK, (min >> MANTISSA) - 1);
    }

    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
                
//        return (float) pow(1 - exp((float)(-numHash * counts.popCount()) / size), numHash);
        return (float) pow((double)(counts.popCount()) / (double)(size), numHash);
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
        this.counts.empty();
    }
    
    public void destroy() {
        this.counts.destroy();
        this.counts = null;
    }
    
    public boolean equivalent(CountingBloomFilter bf) {
        return this.size == bf.size &&
                this.numHash == bf.numHash &&
                BufferComparator.equivalentByteBuffers(counts, bf.counts);
    }
}
