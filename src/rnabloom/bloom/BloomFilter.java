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
import rnabloom.bloom.buffer.LargeBitBuffer;
import rnabloom.bloom.buffer.UnsafeBitBuffer;
import rnabloom.bloom.buffer.AbstractLargeBitBuffer;
import static java.lang.Math.pow;
import static java.lang.Math.exp;

/**
 *
 * @author kmnip
 */
public class BloomFilter implements BloomFilterInterface {    
    protected AbstractLargeBitBuffer bitArray;
    protected final int numHash;
    protected final long size;
    protected final HashFunction hashFunction;
        
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
    
    public static BloomFilter restore(File desc, File bits, HashFunction hashFunction) throws FileNotFoundException, IOException {
        int numHash = -1;
        long size = -1;
        
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
        
        BloomFilter bf = new BloomFilter(size, numHash, hashFunction);
        FileInputStream fin = new FileInputStream(bits);
        bf.bitArray.read(fin);
        fin.close();
        
        return bf;
    }
    
    public void save(File desc, File bits) throws IOException {
        FileWriter writer = new FileWriter(desc);
        
        writer.write(LABEL_SIZE + LABEL_SEPARATOR + this.size + "\n" +
                    LABEL_NUM_HASH + LABEL_SEPARATOR + this.numHash + "\n");
        writer.close();
        
        FileOutputStream out = new FileOutputStream(bits);
        this.bitArray.write(out);
        out.close();
    }
        
    @Override
    public void add(final String key) {
        add(hashFunction.getHashValues(key, numHash));
    }
    
    public void add(final long[] hashVals){
        for (int h=0; h<numHash; ++h) {
            bitArray.set(hashVals[h] % size);
        }
    }

    @Override
    public boolean lookup(final String key) {        
        return lookup(hashFunction.getHashValues(key, numHash));
    }

    public boolean lookup(final long[] hashVals) {
        for (int h=0; h<numHash; ++h) {
            if (!bitArray.get(hashVals[h] % size)) {
                return false;
            }
        }
        
        return true;
    }
    
    @Override
    public float getFPR() {
        /* (1 - e(-kn/m))^k
        k = num hash
        m = size
        n = pop count
        */
        
        return (float) pow(1 - exp(-numHash * bitArray.popCount() / size), numHash);
    }
    
    public void destroy() {
        this.bitArray.destroy();
    }
}
