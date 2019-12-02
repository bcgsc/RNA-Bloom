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
import rnabloom.bloom.hash.CanonicalLeftVariantsNTHashIterator;
import rnabloom.bloom.hash.CanonicalPredecessorsNTHashIterator;
import rnabloom.bloom.hash.CanonicalRightVariantsNTHashIterator;
import rnabloom.bloom.hash.CanonicalSuccessorsNTHashIterator;
import static rnabloom.bloom.hash.HashFunction.combineHashValues;
import static rnabloom.util.SeqUtils.getAltNucleotides;
import static rnabloom.util.SeqUtils.shiftLeft;
import static rnabloom.util.SeqUtils.shiftRight;
import static rnabloom.util.SeqUtils.stringToBytes;
import static rnabloom.util.SeqUtils.getAltNucleotides;

/**
 *
 * @author Ka Ming Nip
 */
public class CanonicalKmer extends Kmer {
    
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
    public long getReverseComplementHash() {
        return getHash();
    }
    
    @Override
    public long getKmerPairHashValue(Kmer rightPartner) {
        if (rightPartner instanceof CanonicalKmer) {
            return getKmerPairHashValue((CanonicalKmer) rightPartner);
        }
        
        return combineHashValues(this.fHashVal, rightPartner.fHashVal);
    }
    
    public long getKmerPairHashValue(CanonicalKmer rightPartner) {
        return Math.min(combineHashValues(fHashVal, rightPartner.fHashVal),
                            combineHashValues(rightPartner.rHashVal, rHashVal));
    }
    
//    @Override
//    public int hashCode(){
//        return (int) getHash(); 
//    }
        
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
        itr.start(fHashVal, rHashVal, bytes[k-1]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
            if (graph.contains(hVals)) {
                return true;
            }
        }
        
        return false;
    }
    
    @Override
    public boolean hasSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, bytes[0]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
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
        itr.start(fHashVal, rHashVal, bytes[k-1]);
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
    
    @Override
    public boolean hasAtLeastXSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, int x) {
        int result = 0;

        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, bytes[0]);
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
    
    @Override
    public int getNumPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        int result = 0;
        
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, bytes[k-1]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
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
        itr.start(fHashVal, rHashVal, bytes[0]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
            if (graph.contains(hVals)) {
                ++result;
            }
        }
        
        return result;
    }
    
    @Override
    public ArrayDeque<Kmer> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        return getPredecessors(k, numHash, graph, 1);
    }
    
    @Override
    public ArrayDeque<Kmer> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, float minKmerCov) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
        
        getPredecessors(k, numHash, graph, result, minKmerCov);
        
        return result;
    }
    
    @Override
    public void getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, ArrayDeque<Kmer> result, float minKmerCov) {
  
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, bytes[k-1]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
            float myCount = graph.getCount(hVals);
            if (myCount >= minKmerCov) {
                byte[] myBytes = shiftRight(this.bytes, k);
                myBytes[0] = itr.currentChar();
                result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
            }
        }
    }
    
    @Override
    public ArrayDeque<Kmer> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        return getSuccessors(k, numHash, graph, 1);
    }

    @Override
    public ArrayDeque<Kmer> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, float minKmerCov) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
        
        getSuccessors(k, numHash, graph, result, minKmerCov);
        
        return result;
    }
    
    @Override
    public void getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, ArrayDeque<Kmer> result, float minKmerCov) {
        
        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, bytes[0]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
            float myCount = graph.getCount(hVals);
            if (myCount >= minKmerCov) {
                byte[] myBytes = shiftLeft(this.bytes, k);
                myBytes[k-1] = itr.currentChar();
                result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
            }
        }
    }
    
    @Override
    public ArrayDeque<Kmer> getPredecessors(int k, int numHash, BloomFilterDeBruijnGraph graph, BloomFilter bf) {        
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
                       
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, bytes[k-1]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
            if (bf.lookup(hVals)) {
                float myCount = graph.getCount(hVals);
                if (myCount > 0) {
                    byte[] myBytes = shiftRight(this.bytes, k);
                    myBytes[0] = itr.currentChar();
                    result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
                }
            }
        }
        
        return result;
    }
            
    @Override
    public ArrayDeque<Kmer> getSuccessors(int k, int numHash, BloomFilterDeBruijnGraph graph, BloomFilter bf) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);

        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, bytes[0]);
        long[] hVals = itr.hVals;
        
        while (itr.hasNext()) {
            itr.next();
            if (bf.lookup(hVals)) {
                float myCount = graph.getCount(hVals);
                if (myCount > 0) {
                    byte[] myBytes = shiftLeft(this.bytes, k);
                    myBytes[k-1] = itr.currentChar();
                    result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
                }
            }
        }
        
        return result;
    }
    
    @Override
    public ArrayDeque<Kmer> getLeftVariants(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        return getLeftVariants(k, numHash, graph, 1);
    }
    
    @Override
    public ArrayDeque<Kmer> getLeftVariants(int k, int numHash, BloomFilterDeBruijnGraph graph, float minKmerCov) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
        
        CanonicalLeftVariantsNTHashIterator itr = new CanonicalLeftVariantsNTHashIterator(k, numHash);
        byte charOut = bytes[0];
        itr.start(fHashVal, rHashVal, charOut);
        long[] hVals = itr.hVals;
        
        for (byte charIn : getAltNucleotides(charOut)) {
            itr.next(charIn);
            float myCount = graph.getCount(hVals);
            if (myCount >= minKmerCov) {
                byte[] myBytes = Arrays.copyOf(bytes, k);
                myBytes[0] = charIn;
                result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
            }
        }
                
        return result;
    }
    
    
    @Override
    public ArrayDeque<Kmer> getRightVariants(int k, int numHash, BloomFilterDeBruijnGraph graph) {
        return getRightVariants(k, numHash, graph, 1);
    }
    
    @Override
    public ArrayDeque<Kmer> getRightVariants(int k, int numHash, BloomFilterDeBruijnGraph graph, float minKmerCov) {
        ArrayDeque<Kmer> result = new ArrayDeque<>(4);
        
        CanonicalRightVariantsNTHashIterator itr = new CanonicalRightVariantsNTHashIterator(k, numHash);
        int kMinus1 = k-1;
        byte charOut = bytes[kMinus1];
        itr.start(fHashVal, rHashVal, charOut);
        long[] hVals = itr.hVals;
        
        for (byte charIn : getAltNucleotides(charOut)) {
            itr.next(charIn);
            float myCount = graph.getCount(hVals);
            if (myCount >= minKmerCov) {
                byte[] myBytes = Arrays.copyOf(bytes, k);
                myBytes[kMinus1] = charIn;
                result.add(new CanonicalKmer(myBytes, myCount, itr.fHashVal, itr.rHashVal));
            }
        }
                
        return result;
    }
    
    @Override
    public boolean hasDepthRight(int k, int numHash, BloomFilterDeBruijnGraph graph, int depth) {
        ArrayDeque<CanonicalSuccessorsNTHashIterator> stack = new ArrayDeque<>(depth);
        CanonicalSuccessorsNTHashIterator itr = new CanonicalSuccessorsNTHashIterator(k, numHash);
        itr.start(fHashVal, rHashVal, bytes[0]);
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
                extension[extensionLength-1] = itr.currentChar();
                
                CanonicalSuccessorsNTHashIterator nextItr = new CanonicalSuccessorsNTHashIterator(k, numHash);
                
                if (extensionLength < k) {
                    nextItr.start(itr.fHashVal, itr.rHashVal, bytes[extensionLength]);
                }
                else {
                    nextItr.start(itr.fHashVal, itr.rHashVal, extension[extensionLength-k]);
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
    
    @Override
    public boolean hasDepthLeft(int k, int numHash, BloomFilterDeBruijnGraph graph, int depth) {
        ArrayDeque<CanonicalPredecessorsNTHashIterator> stack = new ArrayDeque<>(depth);
        CanonicalPredecessorsNTHashIterator itr = new CanonicalPredecessorsNTHashIterator(k, numHash);
        
        itr.start(fHashVal, rHashVal, bytes[k-1]);
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
                extension[extensionLength-1] = itr.currentChar();
                
                CanonicalPredecessorsNTHashIterator nextItr = new CanonicalPredecessorsNTHashIterator(k, numHash);
                
                if (extensionLength < k) {
                    nextItr.start(itr.fHashVal, itr.rHashVal, bytes[k-1-extensionLength]);
                }
                else {
                    nextItr.start(itr.fHashVal, itr.rHashVal, extension[extensionLength-k]);
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
