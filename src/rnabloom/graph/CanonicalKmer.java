/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.graph;

import java.util.ArrayDeque;
import java.util.Arrays;
import rnabloom.bloom.BloomFilter;
import rnabloom.bloom.hash.CanonicalLeftVariantsNTHashIterator;
import rnabloom.bloom.hash.CanonicalPredecessorsNTHashIterator;
import rnabloom.bloom.hash.CanonicalRightVariantsNTHashIterator;
import rnabloom.bloom.hash.CanonicalSuccessorsNTHashIterator;
import static rnabloom.bloom.hash.HashFunction2.combineHashValues;
import static rnabloom.bloom.hash.NTHash.NTM64;
import static rnabloom.util.SeqUtils.NUCLEOTIDES;
import static rnabloom.util.SeqUtils.getAltNucleotides;
import static rnabloom.util.SeqUtils.shiftLeft;
import static rnabloom.util.SeqUtils.shiftRight;
import static rnabloom.util.SeqUtils.stringToBytes;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalKmer extends Kmer2 {
    
    protected long rHashVal;
    
    public CanonicalKmer (String seq, int k, float count, long fHashVal, long rHashVal) {
        this(stringToBytes(seq, k), count, fHashVal, rHashVal);
    }
    
    public CanonicalKmer (byte[] bytes, float count, long fHashVal, long rHashVal) {
        super(bytes, count, fHashVal); // copy bytes
        this.rHashVal = rHashVal;
    }
    
    @Override
    public long getHash() {
        return Math.min(fHashVal, rHashVal);
    }
    
    @Override
    public long[] getHashValues(int len, int numHash) {
        long[] hVals = new long[numHash];
        NTM64(getHash(), hVals, len, numHash);
        return hVals;
    }
    
    @Override
    public long[] getKmerPairHashValues(int len, int numHash, Kmer2 rightPartner) {
        long[] hVals = new long[numHash];
        if (rightPartner instanceof CanonicalKmer) {
            NTM64(Math.min(combineHashValues(fHashVal, rightPartner.fHashVal), combineHashValues(((CanonicalKmer) rightPartner).rHashVal, rHashVal)), hVals, len, numHash);
        }
        else {
            NTM64(combineHashValues(this.fHashVal, rightPartner.fHashVal), hVals, len, numHash);
        }
        return hVals;
    }
    
    public long[] getKmerPairHashValues(int len, int numHash, CanonicalKmer rightPartner) {
        long[] hVals = new long[numHash];
        NTM64(Math.min(combineHashValues(fHashVal, rightPartner.fHashVal), combineHashValues(rightPartner.rHashVal, rHashVal)), hVals, len, numHash);
        return hVals;
    }
    
    public long[] getFHashValues(int len, int numHash) {
        long[] hVals = new long[numHash];
        NTM64(fHashVal, hVals, len, numHash);
        return hVals;
    }
    
    public long[] getRHashValues(int len, int numHash) {
        long[] hVals = new long[numHash];
        NTM64(rHashVal, hVals, len, numHash);
        return hVals;
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
        return Arrays.equals(bytes, ((CanonicalKmer)obj).bytes);
    }
    
    public long getFHash() {
        return fHashVal;
    }
    
    public long getRHash() {
        return rHashVal;
    }
    
    @Override
    public boolean hasPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                return true;
            }
        }
        
        return false;
    }
    
    @Override
    public boolean hasSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                return true;
            }
        }

        return false;
    }
    
    @Override
    public boolean hasAtLeastXPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, int x) {
        int result = 0;
        
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[k-1]);
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
    
    @Override
    public boolean hasAtLeastXSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, int x) {
        int result = 0;

        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[0]);
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
    
    @Override
    public int getNumPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        int result = 0;
        
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                ++result;
            }
        }
        
        return result;
    }
    
    @Override
    public int getNumSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        int result = 0;

        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (graph.contains(hVals)) {
                ++result;
            }
        }
        
        return result;
    }
    
    @Override
    public ArrayDeque<Kmer2> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
  
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        float myCount;
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            myCount = graph.getCount(hVals);
            if (myCount > 0) {
                byte[] myBytes = shiftRight(this.bytes, k);
                myBytes[0] = (byte) c;
                result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
            }
        }
        
        return result;
    }
    
    @Override
    public ArrayDeque<Kmer2> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {        
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
        
        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        float myCount;
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            myCount = graph.getCount(hVals);
            if (myCount > 0) {
                byte[] myBytes = shiftLeft(this.bytes, k);
                myBytes[k-1] = (byte) c;
                result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
            }
        }
        
        return result;
    }
    
    @Override
    public ArrayDeque<Kmer2> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, BloomFilter bf) {        
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
                       
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        float myCount;
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (bf.lookup(hVals)) {
                myCount = graph.getCount(hVals);
                if (myCount > 0) {
                    byte[] myBytes = shiftRight(this.bytes, k);
                    myBytes[0] = (byte) c;
                    result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
                }
            }
        }
        
        return result;
    }
            
    @Override
    public ArrayDeque<Kmer2> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, BloomFilter bf) {
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);

        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        float myCount;
        for (char c : NUCLEOTIDES) {
            itr.next(c);
            if (bf.lookup(hVals)) {
                myCount = graph.getCount(hVals);
                if (myCount > 0) {
                    byte[] myBytes = shiftLeft(this.bytes, k);
                    myBytes[k-1] = (byte) c;
                    result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
                }
            }
        }
        
        return result;
    }
    
    @Override
    public ArrayDeque<Kmer2> getLeftVariants(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
        
        CanonicalLeftVariantsNTHashIterator itr = new CanonicalLeftVariantsNTHashIterator(k, numHash);
        char charOut = (char) bytes[0];
        itr.start(fHashVal, rHashVal, charOut);
        long[] hVals = itr.hVals;
        
        byte[] myBytes;
        float myCount;
        for (char charIn : getAltNucleotides(charOut)) {
            itr.next(charIn);
            myCount = graph.getCount(hVals);
            if (myCount > 0) {
                myBytes = Arrays.copyOf(bytes, k);
                myBytes[0] = (byte) charIn;
                result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
            }
        }
                
        return result;
    }
    
    @Override
    public ArrayDeque<Kmer2> getRightVariants(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        ArrayDeque<Kmer2> result = new ArrayDeque<>(4);
        
        CanonicalRightVariantsNTHashIterator itr = new CanonicalRightVariantsNTHashIterator(k, numHash);
        char charOut = (char) bytes[k-1];
        itr.start(fHashVal, rHashVal, charOut);
        long[] hVals = itr.hVals;
        
        byte[] myBytes;
        float myCount;
        for (char charIn : getAltNucleotides(charOut)) {
            itr.next(charIn);
            myCount = graph.getCount(hVals);
            if (myCount > 0) {
                myBytes = Arrays.copyOf(bytes, k);
                myBytes[k-1] = (byte) charIn;
                result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
            }
        }
                
        return result;
    }
}
