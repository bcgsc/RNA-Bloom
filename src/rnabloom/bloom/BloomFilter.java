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
package rnabloom.bloom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import rnabloom.bloom.buffer.LargeBitBuffer;
import rnabloom.bloom.buffer.UnsafeBitBuffer;
import rnabloom.bloom.buffer.AbstractLargeBitBuffer;
import static java.lang.Math.pow;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import rnabloom.bloom.buffer.BufferComparator;
import rnabloom.bloom.hash.HashFunction;

/**
 *
 * @author Ka Ming Nip
 */
public class BloomFilter implements BloomFilterInterface {    
    protected AbstractLargeBitBuffer bitArray;
    protected int numHash;
    protected long size;
    protected HashFunction hashFunction;
    protected long popcount = -1;
        
    public BloomFilter(long size, int numHash, HashFunction hashFunction) {
        
        this.size = size;
        try {
            //System.out.println("unsafe");
            this.bitArray = new UnsafeBitBuffer(size);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArray = new LargeBitBuffer(size);
        }
        this.numHash = numHash;
        this.hashFunction = hashFunction;
    }
    
    private static final String LABEL_SEPARATOR = ":";
    private static final String LABEL_SIZE = "size";
    private static final String LABEL_NUM_HASH = "numhash";
    private static final String LABEL_FPR = "fpr";
    
    public BloomFilter(File desc, File bits, HashFunction hashFunction) throws FileNotFoundException, IOException {
        this(desc, bits, hashFunction, true);
    }
    
    public BloomFilter(File desc, File bits, HashFunction hashFunction, boolean loadBits) throws FileNotFoundException, IOException {
        
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
            this.bitArray = new UnsafeBitBuffer(size);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArray = new LargeBitBuffer(size);
        }
        
        if (loadBits) {
            FileInputStream fin = new FileInputStream(bits);
            this.bitArray.read(fin);
            fin.close();
        }
        
        /**@TODO Assert file size*/
    }
    
    protected long getIndex(long hashVal) {
        // shift right to remove sign bit and modulus the size of buffer
        return (hashVal >>> 1) % size;
    }
    
    public void save(File desc, File bits) throws IOException {
        FileWriter writer = new FileWriter(desc, false);
        
        writer.write(LABEL_SIZE + LABEL_SEPARATOR + this.size + "\n" +
                    LABEL_NUM_HASH + LABEL_SEPARATOR + this.numHash + "\n" +
                    LABEL_FPR + LABEL_SEPARATOR + this.getFPR() + "\n");
        writer.close();
        
        FileOutputStream out = new FileOutputStream(bits, false);
        this.bitArray.write(out);
        out.close();
    }
        
    @Override
    public void add(String key) {
        final long[] hashVals = new long[numHash];
        hashFunction.getHashValues(key, numHash, hashVals);
        add(hashVals);
    }
    
    public void add(final long[] hashVals){
        for (int h=0; h<numHash; ++h) {
            bitArray.set(getIndex(hashVals[h]));
        }
    }
    
    public void add(final long hashVal) {
        add(hashFunction.getHashValues(hashVal, numHash));
    }
    
    public boolean lookupThenAdd(final long[] hashVals) {
        boolean found = true;
        
        for (int h=0; h<numHash; ++h) {
            found = bitArray.getAndSet(getIndex(hashVals[h])) && found;
        }
        
        return found;
    }
    
    public void addCAS(final long[] hashVals) {
        for (int h=0; h<numHash; ++h) {
            bitArray.setCAS(getIndex(hashVals[h]));
        }        
    }

    @Override
    public boolean lookup(String key) {
        final long[] hashVals = new long[numHash];
        hashFunction.getHashValues(key, numHash, hashVals);
        return lookup(hashVals);
    }

    public boolean lookup(final long[] hashVals) {
        for (int h=0; h<numHash; ++h) {
            if (!bitArray.get(getIndex(hashVals[h]))) {
                return false;
            }
        }
        
        return true;
    }
    
    public boolean lookup(final long hashVal) {
        return lookup(hashFunction.getHashValues(hashVal, numHash));
    }
    
    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
        popcount = bitArray.popCount();
        return (float) pow((double)(popcount) / (double)(size), numHash);
    }

    public static long getExpectedSize(long expNumElements, float fpr, int numHash) {
        double r = (double) (-numHash) / log(1 - exp(log(fpr) / (double) numHash));
        return (long) Math.ceil(expNumElements * r);
    }

    public long getPopCount() {
        return popcount;
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
    
//    public long updatePopcount() {
//        popcount = bitArray.popCount();
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
    
    public boolean equivalent(BloomFilter bf) {
        return this.size == bf.size && 
                this.numHash == bf.numHash &&
                BufferComparator.equivalentBitBuffers(bitArray, bf.bitArray);
    }
}
