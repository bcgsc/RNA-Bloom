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
 * @author gengar
 */
public class PairedKeysPartitionedBloomFilter {

    protected AbstractLargeBitBuffer bitArrayLeft;
    protected AbstractLargeBitBuffer bitArrayRight;
    protected AbstractLargeBitBuffer bitArrayPair;
    protected int numHash;
    protected long size;
    protected long partitionSize;
    protected HashFunction hashFunction;
    
    public PairedKeysPartitionedBloomFilter(long size, int numHash, HashFunction hashFunction) {
        this.size = size;
        this.partitionSize = size/3;
        try {
            //System.out.println("unsafe");
            this.bitArrayLeft = new UnsafeBitBuffer(partitionSize);
            this.bitArrayRight = new UnsafeBitBuffer(partitionSize);
            this.bitArrayPair = new UnsafeBitBuffer(partitionSize);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArrayLeft = new LargeBitBuffer(partitionSize);
            this.bitArrayRight = new LargeBitBuffer(partitionSize);
            this.bitArrayPair = new LargeBitBuffer(partitionSize);
        }
        this.numHash = numHash;
        this.hashFunction = hashFunction;
    }
    
    private static final String LABEL_SEPARATOR = ":";
    private static final String LABEL_SIZE = "size";
    private static final String LABEL_PARTITION_SIZE = "partitionSize";
    private static final String LABEL_NUM_HASH = "numhash";
    private static final String LABEL_FPR = "fpr";
    private static final String LABEL_FPR_LEFT = "fprLeft";
    private static final String LABEL_FPR_RIGHT = "fprRight";
    private static final String LABEL_FPR_PAIR = "fprPair";
    
    public PairedKeysPartitionedBloomFilter(File desc, 
            File leftBits, File rightBits, File pairBits, 
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
                    partitionSize = size/3;
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
            this.bitArrayLeft = new UnsafeBitBuffer(partitionSize);
            this.bitArrayRight = new UnsafeBitBuffer(partitionSize);
            this.bitArrayPair = new UnsafeBitBuffer(partitionSize);
        }
        catch(NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            this.bitArrayLeft = new LargeBitBuffer(partitionSize);
            this.bitArrayRight = new LargeBitBuffer(partitionSize);
            this.bitArrayPair = new LargeBitBuffer(partitionSize);
        }
        
        FileInputStream fin = new FileInputStream(leftBits);
        this.bitArrayLeft.read(fin);
        fin.close();
        
        fin = new FileInputStream(rightBits);
        this.bitArrayRight.read(fin);
        fin.close();
        
        fin = new FileInputStream(pairBits);
        this.bitArrayPair.read(fin);
        fin.close();
        
        /**@TODO Assert file size*/
    }

    public int getNumhash() {
        return numHash;
    }
    
    public void save(File desc, File leftBits, File rightBits, File pairBits) throws IOException {
        FileWriter writer = new FileWriter(desc, false);
        
        float leftFPR = this.getLeftFPR();
        float rightFPR = this.getRightFPR();
        float pairFPR = this.getPairFPR();
        
        writer.write(LABEL_SIZE + LABEL_SEPARATOR + this.size + "\n" +
                    LABEL_PARTITION_SIZE + LABEL_SEPARATOR + this.partitionSize + "\n" +
                    LABEL_NUM_HASH + LABEL_SEPARATOR + this.numHash + "\n" +
                    LABEL_FPR + LABEL_SEPARATOR + (leftFPR * rightFPR * pairFPR) + "\n" +
                    LABEL_FPR_LEFT + LABEL_SEPARATOR + leftFPR + "\n" +
                    LABEL_FPR_RIGHT + LABEL_SEPARATOR + rightFPR + "\n" +
                    LABEL_FPR_PAIR + LABEL_SEPARATOR + pairFPR + "\n"
                );
        writer.close();
        
        FileOutputStream out = new FileOutputStream(leftBits, false);
        this.bitArrayLeft.write(out);
        out.close();
        
        out = new FileOutputStream(rightBits, false);
        this.bitArrayRight.write(out);
        out.close();
        
        out = new FileOutputStream(pairBits, false);
        this.bitArrayPair.write(out);
        out.close();
    }

    protected long getIndex(long hashVal) {
        // shift right to remove sign bit and modulus the size of buffer
        return (hashVal >>> 1) % partitionSize;
    }
    
//    public void add(String left, String right) {
//        long[] hashValsLeft = new long[numHash];
//        hashFunction.getHashValues(left, numHash, hashValsLeft);
//        long[] hashValsRight = new long[numHash];
//        hashFunction.getHashValues(right, numHash, hashValsRight);
//        long[] hashValsPair = hashFunction.getHashValues(left, right, numHash);
//        
//        for (int h=0; h<numHash; ++h) {
//            bitArrayLeft.set(getIndex(hashValsLeft[h]));
//            bitArrayRight.set(getIndex(hashValsRight[h]));
//            bitArrayPair.set(getIndex(hashValsPair[h]));
//        }
//    }
    
    public void add(final long hashValLeft, final long hashValRight, final long hashValPair) {        
        add(hashFunction.getHashValues(hashValLeft, numHash),
            hashFunction.getHashValues(hashValRight, numHash),
            hashFunction.getHashValues(hashValPair, numHash));
    }
    
    public void add(final long[] hashValsLeft, final long[] hashValsRight, final long[] hashValsPair) {
        for (int h=0; h<numHash; ++h) {
            bitArrayLeft.set(getIndex(hashValsLeft[h]));
            bitArrayRight.set(getIndex(hashValsRight[h]));
            bitArrayPair.set(getIndex(hashValsPair[h]));
        }
    }
    
//    public void addLeft(final long[] hashVals) {
//        for (int h=0; h<numHash; ++h) {
//            bitArrayLeft.set(getIndex(hashVals[h]));
//        }
//    }
//    
//    public void addRight(final long[] hashVals) {
//        for (int h=0; h<numHash; ++h) {
//            bitArrayRight.set(getIndex(hashVals[h]));
//        }
//    }
//    
//    public void addPair(final long[] hashVals) {
//        for (int h=0; h<numHash; ++h) {
//            bitArrayPair.set(getIndex(hashVals[h]));
//        }
//    }
    
    public boolean lookup(final long hashValLeft, final long hashValRight, final long hashValPair) {
        return lookupLeft(hashValLeft) && lookupRight(hashValRight) && lookupPair(hashValPair);
    }
    
    public boolean lookup(final long[] hashValsLeft, final long[] hashValsRight, final long[] hashValsPair) {
        return lookupLeft(hashValsLeft) && lookupRight(hashValsRight) && lookupPair(hashValsPair);
    }
        
    public boolean lookupLeft(final long[] hashVals) {
        for (int h=0; h<numHash; ++h) {
            if (!bitArrayLeft.get(getIndex(hashVals[h]))) {
                return false;
            }
        }
        
        return true;
    }
    
    public boolean lookupRight(final long[] hashVals) {
        for (int h=0; h<numHash; ++h) {
            if (!bitArrayRight.get(getIndex(hashVals[h]))) {
                return false;
            }
        }
        
        return true;
    }
    
    public boolean lookupPair(final long[] hashVals) {
        for (int h=0; h<numHash; ++h) {
            if (!bitArrayPair.get(getIndex(hashVals[h]))) {
                return false;
            }
        }
        
        return true;
    }

    public boolean lookupLeft(final long hashVal) {
        return lookupLeft(hashFunction.getHashValues(hashVal, numHash));
    }
    
    public boolean lookupRight(final long hashVal) {
        return lookupRight(hashFunction.getHashValues(hashVal, numHash));
    }

    public boolean lookupPair(final long hashVal) {
        return lookupPair(hashFunction.getHashValues(hashVal, numHash));
    }
    
//    public boolean lookupAndAdd(final long[] hashValsLeft, final long[] hashValsRight) {
//        boolean found = true;
//        
//        for (int h=0; h<numHash; ++h) {
//            if (!bitArrayLeft.getAndSet(getIndex(hashValsLeft[h]))) {
//                found = false;
//            }
//        }
//        
//        for (int h=0; h<numHash; ++h) {
//            if (!bitArrayRight.getAndSet(getIndex(hashValsRight[h]))) {
//                found = false;
//            }
//        }
//        
//        long[] hashVals = hashFunction.getHashValues(hashValsLeft, hashValsRight, numHash);
//        
//        for (int h=0; h<numHash; ++h) {
//            if (!bitArrayPair.getAndSet(getIndex(hashVals[h]))) {
//                found = false;
//            }
//        }
//        
//        return found;
//    }
    
    public void destroy() {
        this.bitArrayLeft.destroy();
        this.bitArrayLeft = null;
        
        this.bitArrayRight.destroy();
        this.bitArrayRight = null;
        
        this.bitArrayPair.destroy();
        this.bitArrayPair = null;
    }
    
    public void empty() {
        this.bitArrayLeft.empty();
        this.bitArrayRight.empty();
        this.bitArrayPair.empty();
    }
    
    public boolean equivalent(PairedKeysPartitionedBloomFilter bf) {
        return this.size == bf.size && 
                this.numHash == bf.numHash && 
                BufferComparator.equivalentBitBuffers(this.bitArrayLeft, bf.bitArrayLeft) &&
                BufferComparator.equivalentBitBuffers(this.bitArrayRight, bf.bitArrayRight) &&
                BufferComparator.equivalentBitBuffers(this.bitArrayPair, bf.bitArrayPair);
    }
    
    public float getFPR() {
        return getLeftFPR() * getRightFPR() * getPairFPR();
    }
    
    public float getLeftFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
//        return (float) pow(1 - exp((float)(-numHash * bitArrayLeft.popCount()) / partitionSize), numHash);
        return (float) pow((double)(bitArrayLeft.popCount()) / (double)(size), numHash);
    }
    
    public float getRightFPR() {        
//        return (float) pow(1 - exp((float)(-numHash * bitArrayRight.popCount()) / partitionSize), numHash);
        return (float) pow((double)(bitArrayRight.popCount()) / (double)(size), numHash);
    }
    
    public float getPairFPR() {        
//        return (float) pow(1 - exp((float)(-numHash * bitArrayPair.popCount()) / partitionSize), numHash);
        return (float) pow((double)(bitArrayPair.popCount()) / (double)(size), numHash);
    }
}
