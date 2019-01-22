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
package rnabloom.graph;

import java.util.ArrayDeque;
import java.util.Arrays;
import rnabloom.bloom.BloomFilter;
import static rnabloom.bloom.hash.HashFunction.combineHashValues;
import rnabloom.bloom.hash.LeftVariantsNTHashIterator;
import rnabloom.bloom.hash.PredecessorsNTHashIterator;
import rnabloom.bloom.hash.SuccessorsNTHashIterator;
import rnabloom.bloom.hash.RightVariantsNTHashIterator;
import static rnabloom.util.SeqUtils.bytesToString;
import static rnabloom.util.SeqUtils.stringToBytes;
import static rnabloom.util.SeqUtils.shiftLeft;
import static rnabloom.util.SeqUtils.shiftRight;
import static rnabloom.util.SeqUtils.getAltNucleotides;

/**
 *
 * @author Ka Ming Nip
 */
public class Kmer {
    
    public byte[] bytes;
    public float count;
    protected long fHashVal;
    
    public Kmer(String seq, int k, float count, long fHashVal) {
        this(stringToBytes(seq, k), count, fHashVal);
    }
    
    public Kmer(byte[] bytes, float count, long fHashVal) {
        this.bytes = bytes; // copy bytes?
        this.count = count;
        this.fHashVal = fHashVal;
    }
    
    public long getHash() {
        return fHashVal;
    }
    
    public long getKmerPairHashValue(Kmer rightPartner) {
        return combineHashValues(this.fHashVal, rightPartner.fHashVal);
    }
                
    public boolean equals(Kmer other) {
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
        return Arrays.equals(bytes, ((Kmer)obj).bytes);
    }
    
    public boolean hasPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
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
        
        while (itr.hasNext()) {
            itr.next();
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
        
        while (itr.hasNext()) {
            itr.next();
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
        
        while (itr.hasNext()) {
            itr.next();
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
        
        while (itr.hasNext()) {
            itr.next();
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
        
        while (itr.hasNext()) {
            itr.next();
            if (graph.contains(hVals)) {
                ++result;
            }
        }
        
        return result;
    }
    
    public ArrayDeque<Kmer> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        return getPredecessors(k, numHash, graph, 1);
    }
    
    public ArrayDeque<Kmer> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, float minKmerCov) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
        
        getPredecessors(k, numHash, graph, result, minKmerCov);
        
        return result;
    }
    
    public void getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, ArrayDeque<Kmer> result, float minKmerCov) {
  
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        float myCount;
        while (itr.hasNext()) {
            itr.next();
            myCount = graph.getCount(hVals);
            if (myCount >= minKmerCov) {
                byte[] myBytes = shiftRight(this.bytes, k);
                myBytes[0] = (byte) itr.currentChar();
                result.add(new Kmer(myBytes, myCount, hVals[0]));
            }
        }
    }
    
    public ArrayDeque<Kmer> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        return getSuccessors(k, numHash, graph, 1);
    }
    
    public ArrayDeque<Kmer> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, float minKmerCov) {  
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
        
        getSuccessors(k, numHash, graph, result, minKmerCov);
        
        return result;
    }
    
    public void getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, ArrayDeque<Kmer> result, float minKmerCov) {
        
        SuccessorsNTHashIterator itr = new SuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        float myCount;
        while (itr.hasNext()) {
            itr.next();
            myCount = graph.getCount(hVals);
            if (myCount >= minKmerCov) {
                byte[] myBytes = shiftLeft(this.bytes, k);
                myBytes[k-1] = (byte) itr.currentChar();
                result.add(new Kmer(myBytes, myCount, hVals[0]));
            }
        }
    }
    
    public ArrayDeque<Kmer> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, BloomFilter bf) {        
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
                       
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[k-1]);
        long[] hVals = itr.hVals;
        
        float myCount;
        while (itr.hasNext()) {
            itr.next();
            if (bf.lookup(hVals)) {
                myCount = graph.getCount(hVals);
                if (myCount > 0) {
                    byte[] myBytes = shiftRight(this.bytes, k);
                    myBytes[0] = (byte) itr.currentChar();
                    result.add(new Kmer(myBytes, myCount, hVals[0]));
                }
            }
        }
        
        return result;
    }
            
    public ArrayDeque<Kmer> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, BloomFilter bf) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);

        SuccessorsNTHashIterator itr = new SuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[0]);
        long[] hVals = itr.hVals;
        
        float myCount;
        while (itr.hasNext()) {
            itr.next();
            if (bf.lookup(hVals)) {
                myCount = graph.getCount(hVals);
                if (myCount > 0) {
                    byte[] myBytes = shiftLeft(this.bytes, k);
                    myBytes[k-1] = (byte) itr.currentChar();
                    result.add(new Kmer(myBytes, myCount, hVals[0]));
                }
            }
        }
        
        return result;
    }
    
    public ArrayDeque<Kmer> getLeftVariants(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        return getLeftVariants(k, numHash, graph, 1);
    }
    
    public ArrayDeque<Kmer> getLeftVariants(int k, int numHash, BloomFilterDeBruijnGraph graph, float minKmerCov) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
        
        LeftVariantsNTHashIterator itr = new LeftVariantsNTHashIterator(k, numHash);
        char charOut = (char) bytes[0];
        itr.start(fHashVal, charOut);
        long[] hVals = itr.hVals;
        
        byte[] myBytes;
        float myCount;
        for (char charIn : getAltNucleotides(charOut)) {
            itr.next(charIn);
            myCount = graph.getCount(hVals);
            if (myCount >= minKmerCov) {
                myBytes = Arrays.copyOf(bytes, k);
                myBytes[0] = (byte) charIn;
                result.add(new Kmer(myBytes, myCount, hVals[0]));
            }
        }
                
        return result;
    }
    
    public ArrayDeque<Kmer> getRightVariants(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        return getRightVariants(k, numHash, graph, 1);
    }
    
    public ArrayDeque<Kmer> getRightVariants(int k, int numHash, BloomFilterDeBruijnGraph graph, float minKmerCov) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
        
        RightVariantsNTHashIterator itr = new RightVariantsNTHashIterator(k, numHash);
        char charOut = (char) bytes[k-1];
        itr.start(fHashVal, charOut);
        long[] hVals = itr.hVals;
        
        byte[] myBytes;
        float myCount;
        for (char charIn : getAltNucleotides(charOut)) {
            itr.next(charIn);
            myCount = graph.getCount(hVals);
            if (myCount >= minKmerCov) {
                myBytes = Arrays.copyOf(bytes, k);
                myBytes[k-1] = (byte) charIn;
                result.add(new Kmer(myBytes, myCount, hVals[0]));
            }
        }
                
        return result;
    }
    
    public boolean hasDepthRight(int k, int numHash, BloomFilterDeBruijnGraph graph, int depth) {
        ArrayDeque<SuccessorsNTHashIterator> stack = new ArrayDeque<>(depth);
        SuccessorsNTHashIterator itr = new SuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, (char) bytes[0]);
        stack.add(itr);
        
        byte[] extension = new byte[depth];
        int extensionLength = 1;
        while (extensionLength > 0) {
            if (extensionLength >= depth) {
                return true;
            }
            
            itr = stack.getLast();
            
            if (itr.hasNext()) {
                itr.next();
                extension[extensionLength-1] = (byte) itr.currentChar();
                
                SuccessorsNTHashIterator nextItr = new SuccessorsNTHashIterator(k, numHash);
                
                if (extensionLength < k) {
                    nextItr.start(itr.hVals[0], (char) bytes[extensionLength]);
                }
                else {
                    nextItr.start(itr.hVals[0], (char) extension[extensionLength-k]);
                }
                
                stack.addLast(nextItr);
                ++extensionLength;
            }
            else {
                stack.pollLast();
                --extensionLength;
            }
        }
        
        return false;
    }
    
    public boolean hasDepthLeft(int k, int numHash, BloomFilterDeBruijnGraph graph, int depth) {
        ArrayDeque<PredecessorsNTHashIterator> stack = new ArrayDeque<>(depth);
        PredecessorsNTHashIterator itr = new PredecessorsNTHashIterator(k, numHash);
        
        itr.start(fHashVal, (char) bytes[k-1]);
        stack.add(itr);
        
        byte[] extension = new byte[depth];
        int extensionLength = 1;
        while (extensionLength > 0) {
            if (extensionLength >= depth) {
                return true;
            }
            
            itr = stack.getLast();
            
            if (itr.hasNext()) {
                itr.next();
                extension[extensionLength-1] = (byte) itr.currentChar();
                
                PredecessorsNTHashIterator nextItr = new PredecessorsNTHashIterator(k, numHash);
                
                if (extensionLength < k) {
                    nextItr.start(itr.hVals[0], (char) bytes[k-1-extensionLength]);
                }
                else {
                    nextItr.start(itr.hVals[0], (char) extension[extensionLength-k]);
                }
                
                stack.addLast(nextItr);
                ++extensionLength;
            }
            else {
                stack.pollLast();
                --extensionLength;
            }
        }
        
        return false;
    }
}
