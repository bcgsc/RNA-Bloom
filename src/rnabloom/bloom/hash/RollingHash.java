/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;

import static rnabloom.util.SeqUtils.kmerize;

/**
 *
 * @author kmnip
 */
public class RollingHash {
    
    private static final long SEED_A = 0x3c8bfbb395c60474L;
    private static final long SEED_C = 0x3193c18562a02b4cL;
    private static final long SEED_G = 0x20323ed082572324L;
    private static final long SEED_T = 0x295549f54be24456L;
    
    private static final long SEED_MULTIPLIER = 0x90b45d39fb6da1faL;
    private static final long SEED_SHIFTER = 27;
    
    private static final long[] ROLL_A = getRollArray(SEED_A);
    private static final long[] ROLL_C = getRollArray(SEED_C);
    private static final long[] ROLL_G = getRollArray(SEED_G);
    private static final long[] ROLL_T = getRollArray(SEED_T);
    private static final int ROLL_SIZE = Long.SIZE;
    
    private final int k;
    
    public RollingHash(int k) {
        this.k = k;
    }
        
    private static long[] getRollArray(final long seed) {
        long[] arr = new long[ROLL_SIZE];
        
        for (int i=0; i<ROLL_SIZE; ++i) {
            arr[i] = Long.rotateRight(seed, i);
        }
        
        return arr;
    }
    
    private static long getNucleotideHash(char c) {
        switch (c) {
            case 'A':
                return SEED_A;
            case 'C':
                return SEED_C;
            case 'G':
                return SEED_G;
            case 'T':
                return SEED_T;
        }
        return 0;
    }
    
    private static long getNucleotideHash(char c, int pos) {
        switch (c) {
            case 'A':
                return ROLL_A[pos % 64];
            case 'C':
                return ROLL_C[pos % 64];
            case 'G':
                return ROLL_G[pos % 64];
            case 'T':
                return ROLL_T[pos % 64];
        }
        return 0;
    }

    /**
     * Return the hash value of kmer 
     * @param kmer
     * @return 
     */
    public long getHashValue(final CharSequence kmer) {
        long hashVal = 0L;
        
        for (int i=0; i<k; ++i) {
            hashVal ^= getNucleotideHash(kmer.charAt(i), i);
        }
        
        return hashVal;
    }
    
    /**
     * Return the hash value of kmer at offset in seq
     * @param seq
     * @param offset    position of first base of kmer
     * @return 
     */
    public long getHashValue(final CharSequence seq, int offset) {
        long hashVal = 0L;
        
        for (int i=0; i<k; ++i) {
            hashVal ^= getNucleotideHash(seq.charAt(offset+i), i);
        }
        
        return hashVal;
    }
    
    public long getHashValueSeeded(final CharSequence kmer, int seed) {
        long hashVal = getHashValue(kmer);
        
        hashVal *= seed ^ k * SEED_MULTIPLIER;
        hashVal ^= hashVal >> SEED_SHIFTER;
        
        return hashVal;
    }
    
    public long getRightHashValue(long hashVal, char remove, char add) {
        return Long.rotateLeft(hashVal ^ getNucleotideHash(remove), 1) ^ getNucleotideHash(add, k-1);
    }
    
    public long getLeftHashValue(long hashVal, char remove, char add) {
        return getNucleotideHash(add) ^ Long.rotateRight(hashVal ^ getNucleotideHash(remove, k-1), 1);
    }
    
    /**
     * Return the hash values of all sliding kmers
     * @param seq
     * @return 
     */
    public long[] getHashValues(final CharSequence seq) {
        int numKmers = seq.length() - k + 1;
        
        long[] out = new long[numKmers];
        
        long hashVal = getHashValue(seq);
        out[0] = hashVal;
        
        for (int i=1; i<numKmers; ++i) {
            hashVal = getRightHashValue(hashVal, seq.charAt(i-1), seq.charAt(i+k-1));
            out[i] = hashVal;
        }
        
        return out;
    }
    
    public static void main(String[] args) {
        String seq = "ATACGTCGATCGTCGTCGTAGTAGATGCTAGCTGACTGATGCGGATCGACTGACTCATCGTACGTAGATCGACTCACACAACCACACTCGATCGTAGCTGACTGACTGATCGTAGCTAGCTAGCTGATAGCTAGCTGAGTACGTAGCTGACTGATCGATCGTAGCTGACTGACTGTACGCGGC";
        int k = 25;
        //int h = 3;
        
        RollingHash rh = new RollingHash(k);
        long[] hashes = rh.getHashValues(seq);
        
        int i = 0;
        String[] kmers = kmerize(seq, k);
        for (String kmer : kmers) {
            long h = rh.getHashValue(kmer);
            long h2 = hashes[i++];
            
            if (h != h2) {
                System.out.println("uhoh");
                return;
            }
        }
    }
}
