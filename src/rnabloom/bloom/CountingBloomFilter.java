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
import rnabloom.bloom.hash.HashFunction;
import rnabloom.bloom.buffer.UnsafeByteBuffer;
import rnabloom.bloom.buffer.AbstractLargeByteBuffer;
import rnabloom.bloom.buffer.LargeByteBuffer;
import static java.lang.Math.exp;
import static java.lang.Math.pow;
import static java.lang.Math.random;
import static java.lang.Math.scalb;

/**
 *
 * @author kmnip
 */
public class CountingBloomFilter implements CountingBloomFilterInterface {
    protected AbstractLargeByteBuffer counts;
    protected int numHash;
    protected long size;
    protected HashFunction hashFunction;
        
    private static final byte MANTISSA = 3;
    private static final byte MANTI_MASK = 0xFF >> (8 - MANTISSA);
    private static final byte ADD_MASK = 0x80 >> (7 - MANTISSA);
    
    public CountingBloomFilter(long size, int numHash, HashFunction hashFunction) {
        this.size = size;
        try {
            System.out.println("unsafe");
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
        FileInputStream fin = new FileInputStream(bytes);
        this.counts.read(fin);
        fin.close();
        
        /**@TODO Assert file size*/
    }
    
    public void save(File desc, File bytes) throws IOException {
        FileWriter writer = new FileWriter(desc);
        
        writer.write(LABEL_SIZE + LABEL_SEPARATOR + this.size + "\n" +
                    LABEL_NUM_HASH + LABEL_SEPARATOR + this.numHash + "\n");
        writer.close();
        
        FileOutputStream out = new FileOutputStream(bytes);
        this.counts.write(out);
        out.close();
    }
        
    @Override
    public void increment(String key) {
        increment(hashFunction.getHashValues(key, numHash));
    }
    
    public synchronized void increment(long[] hashVals) {
        // find the smallest count at all hash positions
        byte min = counts.get(hashVals[0] % size);
        byte c;
        int h;
        for (h=1; h<numHash; ++h) {
            c = counts.get(hashVals[h] % size);
            if (c < min) {
                min = c;
            }
            if (min == 0) {
                break;
            }
        }
        
        // increment the smallest count
        byte updated = min;
        if (min <= MANTI_MASK) {
            ++updated;
        }
        else {
            int shiftVal = (1 << ((updated >> MANTISSA) - 1));
            if ((int) (random() * Integer.MAX_VALUE) % shiftVal == 0) {
                ++updated;
            }
        }
        
        // update the smallest count only
        long index;
        for (h=0; h<numHash; ++h) {
            index = hashVals[h] % size;

            if (counts.get(index) == min) {
                counts.set(index, updated);
            }
        }
    }

    @Override
    public float getCount(String key) {
        return getCount(hashFunction.getHashValues(key, numHash));
    }
    
    public float getCount(long[] hashVals) {
        // find the smallest count
        byte min = counts.get(hashVals[0] % size);
        byte c;
        for (int h=1; h<numHash; ++h) {
            c = counts.get(hashVals[h] % size);
            if (c < min) {
                min = c;
            }
            if (min == 0) {
                break;
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
                
        return (float) pow(1 - exp(-numHash * counts.popCount() / size), numHash);
    }
 
    public void destroy() {
        this.counts.destroy();
    }
}
