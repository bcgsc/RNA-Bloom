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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.pow;
import rnabloom.bloom.buffer.AbstractLargeBitBuffer;
import rnabloom.bloom.buffer.BufferComparator;
import rnabloom.bloom.buffer.LargeBitBuffer;
import rnabloom.bloom.buffer.UnsafeBitBuffer;
import rnabloom.bloom.hash.HashFunction;

/**
 *
 * @author Ka Ming Nip
 */
public class PairedKeysBloomFilter {
    
    protected AbstractLargeBitBuffer bitArrayPair;
    protected int numHash;
    protected long size;
    protected HashFunction hashFunction;
    protected long popcount = -1;
    
    public PairedKeysBloomFilter(long size, int numHash, HashFunction hashFunction) {
        this.size = size;
        try {
            this.bitArrayPair = new UnsafeBitBuffer(size);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArrayPair = new LargeBitBuffer(size);
        }
        this.numHash = numHash;
        this.hashFunction = hashFunction;
    }
    
    private static final String LABEL_SEPARATOR = ":";
    private static final String LABEL_SIZE = "size";
    private static final String LABEL_NUM_HASH = "numhash";
    private static final String LABEL_FPR = "fpr";
    
    public PairedKeysBloomFilter(File desc, File pairBits, 
            HashFunction hashFunction) throws FileNotFoundException, IOException {
        
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
            this.bitArrayPair = new UnsafeBitBuffer(size);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArrayPair = new LargeBitBuffer(size);
        }
        
        FileInputStream fin = new FileInputStream(pairBits);
        this.bitArrayPair.read(fin);
        fin.close();
        
        /**@TODO Assert file size*/
    }

    public int getNumhash() {
        return numHash;
    }
    
    public void save(File desc, File bits) throws IOException {
        FileWriter writer = new FileWriter(desc, false);
        
        writer.write(LABEL_SIZE + LABEL_SEPARATOR + this.size + "\n" +
                    LABEL_NUM_HASH + LABEL_SEPARATOR + this.numHash + "\n" +
                    LABEL_FPR + LABEL_SEPARATOR + this.getFPR() + "\n"
                );
        writer.close();
        
        FileOutputStream out = new FileOutputStream(bits, false);
        this.bitArrayPair.write(out);
        out.close();
    }

    protected static long getIndex(long hashVal, long size) {
        // shift right to remove sign bit and modulus the size of buffer
        return (hashVal >>> 1) % size;
    }
    
//    public void add(String left, String right) {
//        long[] hashValsPair = hashFunction.getHashValues(left, right, numHash);
//        
//        for (int h=0; h<numHash; ++h) {
//            bitArrayPair.set(getIndex(hashValsPair[h]));
//        }
//    }
    
    public void add(final long hashValPair) {
        add(hashFunction.getHashValues(hashValPair, numHash));
    }
    
    public void add(final long[] hashValsPair) {
        for (int h=0; h<numHash; ++h) {
            bitArrayPair.set(getIndex(hashValsPair[h], size));
        }
    }

    public boolean lookup(final long hashValsPair) {
        return lookup(hashFunction.getHashValues(hashValsPair, numHash));
    }
    
    public boolean lookup(final long[] hashValsPair) {
        for (int h=0; h<numHash; ++h) {
            if (!bitArrayPair.get(getIndex(hashValsPair[h], size))) {
                return false;
            }
        }
        
        return true;
    }
    
    public boolean lookupThenAdd(final long hashValsPair) {
        return lookupThenAdd(hashFunction.getHashValues(hashValsPair, numHash));
    }
    
    public boolean lookupThenAdd(final long[] hashVals) {
        boolean foundAll = true;
        
        for (int h=0; h<numHash; ++h) {
            foundAll &= bitArrayPair.getAndSet(getIndex(hashVals[h], size));
        }
        
        return foundAll;
    }

//    public boolean lookup(final long hashValsLeft, final long hashValsRight) {
//        return lookup(hashFunction.getHashValues(hashValsLeft, numHash), hashFunction.getHashValues(hashValsRight, numHash));
//    }    
//    
//    public boolean lookup(final long[] hashValsLeft, final long[] hashValsRight) {
//        boolean found = true;
//        
//        long[] hashVals = hashFunction.getHashValues(hashValsLeft, hashValsRight, numHash);
//        
//        for (int h=0; h<numHash; ++h) {
//            if (!bitArrayPair.get(getIndex(hashVals[h]))) {
//                found = false;
//            }
//        }
//        
//        return found;
//    }
        
    public void destroy() {
        if (this.bitArrayPair != null) {
            this.bitArrayPair.destroy();
            this.bitArrayPair = null;
        }
    }
    
    public void empty() {
        if (this.bitArrayPair != null) {
            this.bitArrayPair.empty();
        }
    }
    
    public boolean equivalent(PairedKeysPartitionedBloomFilter bf) {
        return this.size == bf.size && 
                this.numHash == bf.numHash && 
                BufferComparator.equivalentBitBuffers(this.bitArrayPair, bf.bitArrayPair);
    }
    
    public float getFPR() {      
        popcount = bitArrayPair.popCount();
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
}
