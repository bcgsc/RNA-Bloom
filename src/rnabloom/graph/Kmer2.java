/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.graph;

import java.util.ArrayDeque;
import java.util.Arrays;
import rnabloom.bloom.BloomFilter;
import static rnabloom.bloom.hash.HashFunction2.combineHashValues;
import rnabloom.bloom.hash.LeftVariantsNTHashIterator;
import rnabloom.bloom.hash.PredecessorsNTHashIterator;
import rnabloom.bloom.hash.SuccessorsNTHashIterator;
import rnabloom.bloom.hash.RightVariantsNTHashIterator;
import static rnabloom.bloom.hash.NTHash.NTM64;
import static rnabloom.util.SeqUtils.NUCLEOTIDES;
import static rnabloom.util.SeqUtils.bytesToString;
import static rnabloom.util.SeqUtils.stringToBytes;
import static rnabloom.util.SeqUtils.getAltNucleotides;
import static rnabloom.util.SeqUtils.shiftLeft;
import static rnabloom.util.SeqUtils.shiftRight;

/**
 *
 * @author Ka Ming Nip
 */
public class Kmer2 {
    
    public byte[] bytes;
    public float count;
    protected long fHashVal;
    
    public Kmer2(String seq, int k, float count, long fHashVal) {
        this(stringToBytes(seq, k), count, fHashVal);
    }
    
    public Kmer2(byte[] bytes, float count, long fHashVal) {
        this.bytes = bytes; // copy bytes?
        this.count = count;
        this.fHashVal = fHashVal;
    }
    
    public long getHash() {
        return fHashVal;
    }
    
    public long[] getHashValues(int len, int numHash) {
        long[] hVals = new long[numHash];
        NTM64(getHash(), hVals, len, numHash);
        return hVals;
    }
    
    public long[] getKmerPairHashValues(int len, int numHash, Kmer2 rightPartner) {
        long[] hVals = new long[numHash];
        NTM64(combineHashValues(this.fHashVal, rightPartner.fHashVal), hVals, len, numHash);
        return hVals;
    }
            
    public boolean equals(Kmer2 other) {
        return Arrays.equals(bytes, other.bytes);
    }

    @Override
    public String toString() {
        return bytesToString(bytes, bytes.length);
    }

    @Override
    public int hashCode(){
        return (int) getHash(); 
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        return Arrays.equals(bytes, ((Kmer2)obj).bytes);
    }
    
    public boolean hasPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                return true;
            }
        }
        
        return false;
    }
    
    public boolean hasSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        SuccessorsNTHashIterator itr = new SuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                return true;
            }
        }

        return false;
    }
    
    public boolean hasAtLeastXPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, int x) {
        int result = 0;
        
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                if (++result >= x) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    public boolean hasAtLeastXSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, int x) {
        int result = 0;

        SuccessorsNTHashIterator itr = new SuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                if (++result >= x) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    public int getNumPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        int result = 0;
        
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                ++result;
            }
        }
        
        return result;
    }
    
    public int getNumSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        int result = 0;

        SuccessorsNTHashIterator itr = new SuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                ++result;
            }
        }
        
        return result;
    }
    
    public ArrayDeque<Kmer2> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
  
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        float myCount;
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            myCount = graph.getCount(hVals);
            if (myCount > 0) {
                byte[] myBytes = shiftRight(this.bytes, k);
                myBytes[0] = (byte) c;
                result.add(new Kmer2(myBytes, myCount, hVals[0]));
            }
        }
        
        return result;
    }
    
    public ArrayDeque<Kmer2> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {        
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
        
        SuccessorsNTHashIterator itr = new SuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        float myCount;
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            myCount = graph.getCount(hVals);
            if (myCount > 0) {
                byte[] myBytes = shiftLeft(this.bytes, k);
                myBytes[k-1] = (byte) c;
                result.add(new Kmer2(myBytes, myCount, hVals[0]));
            }
        }
        
        return result;
    }
    
    public ArrayDeque<Kmer2> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, BloomFilter bf) {        
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
                       
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        float myCount;
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (bf.lookup(hVals)) {
                myCount = graph.getCount(hVals);
                if (myCount > 0) {
                    byte[] myBytes = shiftRight(this.bytes, k);
                    myBytes[0] = (byte) c;
                    result.add(new Kmer2(myBytes, myCount, hVals[0]));
                }
            }
        }
        
        return result;
    }
            
    public ArrayDeque<Kmer2> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, BloomFilter bf) {
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);

        SuccessorsNTHashIterator itr = new SuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        float myCount;
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (bf.lookup(hVals)) {
                myCount = graph.getCount(hVals);
                if (myCount > 0) {
                    byte[] myBytes = shiftLeft(this.bytes, k);
                    myBytes[k-1] = (byte) c;
                    result.add(new Kmer2(myBytes, myCount, hVals[0]));
                }
            }
        }
        
        return result;
    }
    
    public ArrayDeque<Kmer2> getLeftVariants(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
        
        LeftVariantsNTHashIterator itr = new LeftVariantsNTHashIterator(k, numHash);
        char charOut = (char) bytes[0];
        itr.start(fHashVal, charOut);
        long[] hVals = itr.hVals;
        
        byte[] myBytes;
        float myCount;
        for (char charIn : getAltNucleotides(charOut)) {
            itr.next(charIn);
            myCount = graph.getCount(hVals);
            if (myCount > 0) {
                myBytes = Arrays.copyOf(bytes, k);
                myBytes[0] = (byte) charIn;
                result.add(new Kmer2(myBytes, myCount, hVals[0]));
            }
        }
                
        return result;
    }
    
    public ArrayDeque<Kmer2> getRightVariants(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
        
        RightVariantsNTHashIterator itr = new RightVariantsNTHashIterator(k, numHash);
        char charOut = (char) bytes[k-1];
        itr.start(fHashVal, charOut);
        long[] hVals = itr.hVals;
        
        byte[] myBytes;
        float myCount;
        for (char charIn : getAltNucleotides(charOut)) {
            itr.next(charIn);
            myCount = graph.getCount(hVals);
            if (myCount > 0) {
                myBytes = Arrays.copyOf(bytes, k);
                myBytes[k-1] = (byte) charIn;
                result.add(new Kmer2(myBytes, myCount, hVals[0]));
            }
        }
                
        return result;
    }
}
