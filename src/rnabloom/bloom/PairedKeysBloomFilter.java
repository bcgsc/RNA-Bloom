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
import static java.lang.Math.pow;
import rnabloom.bloom.buffer.AbstractLargeBitBuffer;
import rnabloom.bloom.buffer.BufferComparator;
import rnabloom.bloom.buffer.LargeBitBuffer;
import rnabloom.bloom.buffer.UnsafeBitBuffer;
import rnabloom.bloom.hash.HashFunction;

/**
 *
 * @author kmnip
 */
public class PairedKeysBloomFilter {
    
    protected AbstractLargeBitBuffer bitArrayPair;
    protected int numHash;
    protected long size;
    protected HashFunction hashFunction;
    
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

    public void save(File desc, File leftBits, File rightBits, File pairBits) throws IOException {
        FileWriter writer = new FileWriter(desc, false);
        
        writer.write(LABEL_SIZE + LABEL_SEPARATOR + this.size + "\n" +
                    LABEL_NUM_HASH + LABEL_SEPARATOR + this.numHash + "\n" +
                    LABEL_FPR + LABEL_SEPARATOR + this.getFPR() + "\n"
                );
        writer.close();
        
        FileOutputStream out = new FileOutputStream(pairBits, false);
        this.bitArrayPair.write(out);
        out.close();
    }

    protected long getIndex(long hashVal) {
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
    
    public void add(final long[] hashValsPair) {
        for (int h=0; h<numHash; ++h) {
            bitArrayPair.set(getIndex(hashValsPair[h]));
        }
    }

    public boolean lookup(final long hashValsPair) {
        return lookup(hashFunction.getHashValues(hashValsPair, numHash));
    }
    
    public boolean lookup(final long[] hashValsPair) {
        for (int h=0; h<numHash; ++h) {
            if (!bitArrayPair.get(getIndex(hashValsPair[h]))) {
                return false;
            }
        }
        
        return true;
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
        this.bitArrayPair.destroy();
        this.bitArrayPair = null;
    }
    
    public void empty() {
        this.bitArrayPair.empty();
    }
    
    public boolean equivalent(PairedKeysPartitionedBloomFilter bf) {
        return this.size == bf.size && 
                this.numHash == bf.numHash && 
                BufferComparator.equivalentBitBuffers(this.bitArrayPair, bf.bitArrayPair);
    }
    
    public float getFPR() {      
//        return (float) pow(1 - exp((float)(-numHash * bitArrayPair.popCount()) / partitionSize), numHash);
        return (float) pow((double)(bitArrayPair.popCount()) / (double)(size), numHash);
    }
    
}
