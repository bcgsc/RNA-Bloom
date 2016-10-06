/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

/**
 *
 * @author kmnip
 */
public class HashFunction {
    
    protected final int seed;
    protected final int k;
    
    protected final int kMod4;
    protected final int kDiv4;
    protected final int maxNumVals;
    protected final long[] hashStartVals;
    
    private static final long MIX = 0xc6a4a7935bd1e995L; 
    private static final int ROLL = 47;
    
    public HashFunction(int seed, int k, int maxNumVals) {
        this.seed = seed;
        this.k = k;
        this.kMod4 = k % 4;
        this.kDiv4 = k / 4;
        this.maxNumVals = maxNumVals;
        this.hashStartVals = new long[maxNumVals];
        for (int i=0; i<maxNumVals; ++i) {
            hashStartVals[i] = ( (seed+i) & 0xffffffffL ) ^ (k*MIX);
        }
    }
        
    public int getSeed() {
        return seed;
    }
    
    public int getK() {
        return k;
    }
    
    /**
     * Generate a 64 bit hash from CharSequence 
     */
    protected long mmh3SingleVal(final long[] charModVals, long h, int offset) {
        for (int i=0; i<kDiv4; ++i) { 
 
            long v = ( ( charModVals[offset++]       ) | 
                       ( charModVals[offset++] << 16 ) | 
                       ( charModVals[offset++] << 32 ) | 
                       ( charModVals[offset++] << 48 ) ); 
 
            v *= MIX; 
            v ^= v >>> ROLL; 
            v *= MIX; 
 
            h ^= v; 
            h *= MIX; 
        } 
 
        switch (kMod4) { 
            case 3: h ^= charModVals[offset+2] << 32;
                    break;
            case 2: h ^= charModVals[offset+1] << 16;
                    break;
            case 1: h ^= charModVals[offset]; 
                    h *= MIX;
                    break;
        } 
 
        h ^= h >>> ROLL; 
        h *= MIX; 
        h ^= h >>> ROLL;
        
        return h >>> 1; // remove the sign bit
    }
    
    
    /**
     * Generate a 64 bit hash from CharSequence 
     */
    protected long mmh3SingleVal(final CharSequence seq, long h, int offset) {
        for (int i=0; i<kDiv4; ++i) { 
 
            long v = ( ( ( seq.charAt( offset++ ) & 0xffffL )       ) | 
                       ( ( seq.charAt( offset++ ) & 0xffffL ) << 16 ) | 
                       ( ( seq.charAt( offset++ ) & 0xffffL ) << 32 ) | 
                       ( ( seq.charAt( offset++ ) & 0xffffL ) << 48 ) ); 
 
            v *= MIX; 
            v ^= v >>> ROLL; 
            v *= MIX; 
 
            h ^= v; 
            h *= MIX; 
        } 
 
        switch (kMod4) { 
            case 3: h ^= ( seq.charAt( offset + 2 ) & 0xffffL ) << 32;
                    break;
            case 2: h ^= ( seq.charAt( offset + 1 ) & 0xffffL ) << 16;
                    break;
            case 1: h ^= ( seq.charAt( offset     ) & 0xffffL ); 
                    h *= MIX;
                    break;
        } 
 
        h ^= h >>> ROLL; 
        h *= MIX; 
        h ^= h >>> ROLL;
        
        return h >>> 1; // remove the sign bit
    }

    public void getHashValues(final CharSequence kmer,
                             final long[] out) {
        
        for (int h=0; h<maxNumVals; ++h) {
            out[h] = mmh3SingleVal(kmer, hashStartVals[h], 0);
        }
    }
    
    public void getHashValues(final CharSequence kmer,
                             final int numHash,
                             final long[] out) {
        
        for (int h=0; h<numHash; ++h) {
            out[h] = mmh3SingleVal(kmer, hashStartVals[h], 0);
        }
    }
    
    public long[] getHashValues(final CharSequence kmer, final int numHash) {
        final long[] hashVals = new long[numHash];
        //murmurhash3_x64_128(kmer.getBytes(), 0, k, seed, numHash, hashVals);
        getHashValues(kmer, numHash, hashVals);         
        
        return hashVals;
    }
        
    public long[] getHashValues(final CharSequence kmer1, final CharSequence kmer2, final int numHash) {
        final long[] hashVals1 = new long[numHash];
        //murmurhash3_x64_128(kmer1.getBytes(), 0, k, seed, numHash, hashVals1);
        getHashValues(kmer1, numHash, hashVals1);
        
        final long[] hashVals2 = new long[numHash];
        //murmurhash3_x64_128(kmer2.getBytes(), 0, k, seed, numHash, hashVals2);
        getHashValues(kmer2, numHash, hashVals2);
        
        final long[] hashVal = new long[numHash];
        for (int i=0; i<numHash; ++i) {
            hashVal[i] = combineHashValues(hashVals1[i], hashVals2[i]);
        }
        
        return hashVal;
    }
    
    public long[][] getHashValues(final CharSequence seq) {
        int seqLength = seq.length();
        int numKmers = seqLength - k + 1;
        
        long[][] out = new long[numKmers][maxNumVals];
        
        long[] charModVals = new long[seqLength];
        for (int i=0; i<seqLength; ++i) {
            charModVals[i] = seq.charAt(i) & 0xffffL;
        }
        
        for (int i=0; i<numKmers; ++i) {
            for (int j=0; j<maxNumVals; ++j) {
                out[i][j] = mmh3SingleVal(charModVals, hashStartVals[j], i);
            }
        }
        
        return out;
    }
    
    public long[][] getHashValuesReverseComplement(final CharSequence seq) {
        int seqLength = seq.length();
        int numKmers = seqLength - k + 1;
        
        long[][] out = new long[numKmers][maxNumVals];
        
        long[] charModVals = new long[seqLength];
        int iLast = seqLength-1;
        for (int i=0; i<seqLength; ++i) {
            switch(seq.charAt(i)) {
                case 'A':
                    charModVals[iLast-i] = 'T' & 0xffffL;
                    break;
                case 'C':
                    charModVals[iLast-i] = 'G' & 0xffffL;
                    break;
                case 'G':
                    charModVals[iLast-i] = 'C' & 0xffffL;
                    break;
                case 'T':
                    charModVals[iLast-i] = 'A' & 0xffffL;
                    break;
                default:
                    charModVals[iLast-i] = 'N' & 0xffffL;
            }
        }
        
        for (int i=0; i<numKmers; ++i) {
            for (int j=0; j<maxNumVals; ++j) {
                out[i][j] = mmh3SingleVal(charModVals, hashStartVals[j], i);
            }
        }
        
        return out;
    }
    
    public static long combineHashValues(long a, long b) {
        // See: http://stackoverflow.com/a/27952689
        a ^= b + 0x9e3779b9 + (a << 6) + (b >> 2);
        return a < 0 ? -a : a;
    }
}
