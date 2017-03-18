/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom.hash;


/**
 *
 * @author Hamid Mohamadi, Genome Sciences Centre, BC Cancer Agency
 */
public class NTHash {
    
    
    // offset for the complement base in the random seeds table
    private static final int cpOff = 0x07;

    // shift for generating multiple hash values
    private static final int multiShift = 27;

    // seed for generating multiple hash values
    private final static long multiSeed = 0x90b45d39fb6da1faL;

    // 64-bit random seeds corresponding to bases and their complements
    private final static long seedA = 0x3c8bfbb395c60474L;
    private final static long seedC = 0x3193c18562a02b4cL;
    private final static long seedG = 0x20323ed082572324L;
    private final static long seedT = 0x295549f54be24456L;
    private final static long seedN = 0x0000000000000000L;

    private final static long[] vecA = {
        0x3c8bfbb395c60474L,0x7917f7672b8c08e8L,0xf22feece571811d0L,0xe45fdd9cae3023a1L,0xc8bfbb395c604743L,0x917f7672b8c08e87L,0x22feece571811d0fL,0x45fdd9cae3023a1eL,
        0x8bfbb395c604743cL,0x17f7672b8c08e879L,0x2feece571811d0f2L,0x5fdd9cae3023a1e4L,0xbfbb395c604743c8L,0x7f7672b8c08e8791L,0xfeece571811d0f22L,0xfdd9cae3023a1e45L,
        0xfbb395c604743c8bL,0xf7672b8c08e87917L,0xeece571811d0f22fL,0xdd9cae3023a1e45fL,0xbb395c604743c8bfL,0x7672b8c08e87917fL,0xece571811d0f22feL,0xd9cae3023a1e45fdL,
        0xb395c604743c8bfbL,0x672b8c08e87917f7L,0xce571811d0f22feeL,0x9cae3023a1e45fddL,0x395c604743c8bfbbL,0x72b8c08e87917f76L,0xe571811d0f22feecL,0xcae3023a1e45fdd9L,
        0x95c604743c8bfbb3L,0x2b8c08e87917f767L,0x571811d0f22feeceL,0xae3023a1e45fdd9cL,0x5c604743c8bfbb39L,0xb8c08e87917f7672L,0x71811d0f22feece5L,0xe3023a1e45fdd9caL,
        0xc604743c8bfbb395L,0x8c08e87917f7672bL,0x1811d0f22feece57L,0x3023a1e45fdd9caeL,0x604743c8bfbb395cL,0xc08e87917f7672b8L,0x811d0f22feece571L,0x023a1e45fdd9cae3L,
        0x04743c8bfbb395c6L,0x08e87917f7672b8cL,0x11d0f22feece5718L,0x23a1e45fdd9cae30L,0x4743c8bfbb395c60L,0x8e87917f7672b8c0L,0x1d0f22feece57181L,0x3a1e45fdd9cae302L,
        0x743c8bfbb395c604L,0xe87917f7672b8c08L,0xd0f22feece571811L,0xa1e45fdd9cae3023L,0x43c8bfbb395c6047L,0x87917f7672b8c08eL,0x0f22feece571811dL,0x1e45fdd9cae3023aL
    };

    private final static long[] vecC = {
        0x3193c18562a02b4cL,0x6327830ac5405698L,0xc64f06158a80ad30L,0x8c9e0c2b15015a61L,0x193c18562a02b4c3L,0x327830ac54056986L,0x64f06158a80ad30cL,0xc9e0c2b15015a618L,
        0x93c18562a02b4c31L,0x27830ac540569863L,0x4f06158a80ad30c6L,0x9e0c2b15015a618cL,0x3c18562a02b4c319L,0x7830ac5405698632L,0xf06158a80ad30c64L,0xe0c2b15015a618c9L,
        0xc18562a02b4c3193L,0x830ac54056986327L,0x06158a80ad30c64fL,0x0c2b15015a618c9eL,0x18562a02b4c3193cL,0x30ac540569863278L,0x6158a80ad30c64f0L,0xc2b15015a618c9e0L,
        0x8562a02b4c3193c1L,0x0ac5405698632783L,0x158a80ad30c64f06L,0x2b15015a618c9e0cL,0x562a02b4c3193c18L,0xac54056986327830L,0x58a80ad30c64f061L,0xb15015a618c9e0c2L,
        0x62a02b4c3193c185L,0xc54056986327830aL,0x8a80ad30c64f0615L,0x15015a618c9e0c2bL,0x2a02b4c3193c1856L,0x54056986327830acL,0xa80ad30c64f06158L,0x5015a618c9e0c2b1L,
        0xa02b4c3193c18562L,0x4056986327830ac5L,0x80ad30c64f06158aL,0x015a618c9e0c2b15L,0x02b4c3193c18562aL,0x056986327830ac54L,0x0ad30c64f06158a8L,0x15a618c9e0c2b150L,
        0x2b4c3193c18562a0L,0x56986327830ac540L,0xad30c64f06158a80L,0x5a618c9e0c2b1501L,0xb4c3193c18562a02L,0x6986327830ac5405L,0xd30c64f06158a80aL,0xa618c9e0c2b15015L,
        0x4c3193c18562a02bL,0x986327830ac54056L,0x30c64f06158a80adL,0x618c9e0c2b15015aL,0xc3193c18562a02b4L,0x86327830ac540569L,0x0c64f06158a80ad3L,0x18c9e0c2b15015a6L
    };

    private final static long[] vecG = {
        0x20323ed082572324L,0x40647da104ae4648L,0x80c8fb42095c8c90L,0x0191f68412b91921L,0x0323ed0825723242L,0x0647da104ae46484L,0x0c8fb42095c8c908L,0x191f68412b919210L,
        0x323ed08257232420L,0x647da104ae464840L,0xc8fb42095c8c9080L,0x91f68412b9192101L,0x23ed082572324203L,0x47da104ae4648406L,0x8fb42095c8c9080cL,0x1f68412b91921019L,
        0x3ed0825723242032L,0x7da104ae46484064L,0xfb42095c8c9080c8L,0xf68412b919210191L,0xed08257232420323L,0xda104ae464840647L,0xb42095c8c9080c8fL,0x68412b919210191fL,
        0xd08257232420323eL,0xa104ae464840647dL,0x42095c8c9080c8fbL,0x8412b919210191f6L,0x08257232420323edL,0x104ae464840647daL,0x2095c8c9080c8fb4L,0x412b919210191f68L,
        0x8257232420323ed0L,0x04ae464840647da1L,0x095c8c9080c8fb42L,0x12b919210191f684L,0x257232420323ed08L,0x4ae464840647da10L,0x95c8c9080c8fb420L,0x2b919210191f6841L,
        0x57232420323ed082L,0xae464840647da104L,0x5c8c9080c8fb4209L,0xb919210191f68412L,0x7232420323ed0825L,0xe464840647da104aL,0xc8c9080c8fb42095L,0x919210191f68412bL,
        0x232420323ed08257L,0x464840647da104aeL,0x8c9080c8fb42095cL,0x19210191f68412b9L,0x32420323ed082572L,0x64840647da104ae4L,0xc9080c8fb42095c8L,0x9210191f68412b91L,
        0x2420323ed0825723L,0x4840647da104ae46L,0x9080c8fb42095c8cL,0x210191f68412b919L,0x420323ed08257232L,0x840647da104ae464L,0x080c8fb42095c8c9L,0x10191f68412b9192L
    };

    private final static long[] vecT = {
        0x295549f54be24456L,0x52aa93ea97c488acL,0xa55527d52f891158L,0x4aaa4faa5f1222b1L,0x95549f54be244562L,0x2aa93ea97c488ac5L,0x55527d52f891158aL,0xaaa4faa5f1222b14L,
        0x5549f54be2445629L,0xaa93ea97c488ac52L,0x5527d52f891158a5L,0xaa4faa5f1222b14aL,0x549f54be24456295L,0xa93ea97c488ac52aL,0x527d52f891158a55L,0xa4faa5f1222b14aaL,
        0x49f54be244562955L,0x93ea97c488ac52aaL,0x27d52f891158a555L,0x4faa5f1222b14aaaL,0x9f54be2445629554L,0x3ea97c488ac52aa9L,0x7d52f891158a5552L,0xfaa5f1222b14aaa4L,
        0xf54be24456295549L,0xea97c488ac52aa93L,0xd52f891158a55527L,0xaa5f1222b14aaa4fL,0x54be24456295549fL,0xa97c488ac52aa93eL,0x52f891158a55527dL,0xa5f1222b14aaa4faL,
        0x4be24456295549f5L,0x97c488ac52aa93eaL,0x2f891158a55527d5L,0x5f1222b14aaa4faaL,0xbe24456295549f54L,0x7c488ac52aa93ea9L,0xf891158a55527d52L,0xf1222b14aaa4faa5L,
        0xe24456295549f54bL,0xc488ac52aa93ea97L,0x891158a55527d52fL,0x1222b14aaa4faa5fL,0x24456295549f54beL,0x488ac52aa93ea97cL,0x91158a55527d52f8L,0x222b14aaa4faa5f1L,
        0x4456295549f54be2L,0x88ac52aa93ea97c4L,0x1158a55527d52f89L,0x22b14aaa4faa5f12L,0x456295549f54be24L,0x8ac52aa93ea97c48L,0x158a55527d52f891L,0x2b14aaa4faa5f122L,
        0x56295549f54be244L,0xac52aa93ea97c488L,0x58a55527d52f8911L,0xb14aaa4faa5f1222L,0x6295549f54be2445L,0xc52aa93ea97c488aL,0x8a55527d52f89115L,0x14aaa4faa5f1222bL
    };

    private final static long[] vecN = {
        seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
        seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
        seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
        seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
        seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
        seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
        seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
        seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN
    };

    private final static long[][] msTab = {
        vecN, vecT, vecN, vecG, vecA, vecN, vecN, vecC, // 0..7
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 8..15
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 16..23
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 24..31
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 32..39
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 40..47
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 48..55
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 56..63
        vecN, vecA, vecN, vecC, vecN, vecN, vecN, vecG, // 64..71
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 72..79
        vecN, vecN, vecN, vecN, vecT, vecN, vecN, vecN, // 80..87
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 88..95
        vecN, vecA, vecN, vecC, vecN, vecN, vecN, vecG, // 96..103
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 104..111
        vecN, vecN, vecN, vecN, vecT, vecN, vecN, vecN, // 112..119
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 120..127
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 128..135
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 136..143
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 144..151
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 152..159
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 160..167
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 168..175
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 176..183
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 184..191
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 192..199
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 200..207
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 208..215
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 216..223
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 224..231
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 232..239
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 240..247
        vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN  // 248..255
    };

    private final static long[] seedTab = {
        seedN, seedT, seedN, seedG, seedA, seedN, seedN, seedC, // 0..7
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 40..47
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 48..55
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
        seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 64..71
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 72..79
        seedN, seedN, seedN, seedN, seedT, seedN, seedN, seedN, // 80..87
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
        seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 96..103
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 104..111
        seedN, seedN, seedN, seedN, seedT, seedN, seedN, seedN, // 112..119
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 120..127
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 128..135
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 136..143
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 144..151
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 152..159
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 160..167
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 168..175
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 176..183
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 184..191
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 192..199
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 200..207
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 208..215
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 216..223
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 224..231
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 232..239
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 240..247
        seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN  // 248..255
    };

//    // rotate "v" to the left by "s" positions
//    public static long rol(final long v, final int s) {
//        return (v << s) | (v >>> (64 - s));
//    }
//
//    // rotate "v" to the right by "s" positions
//    public static long ror(final long v, final int s) {
//        return (v >>> s) | (v << (64 - s));
//    }

    /** 
     * Forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @return          hash value
     */
    private static long getFhval(final CharSequence kmerSeq, final int k) {
        long hVal=0;
        for(int i=0; i<k; ++i)
            hVal ^= Long.rotateLeft(seedTab[kmerSeq.charAt(i)], k-1-i);
        return hVal;
    }

    /**
     * Reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @return          hash value
     */
    private static long getRhval(final CharSequence kmerSeq, final int k) {
        long hVal=0;
        for(int i=0; i<k; ++i)
            hVal ^= Long.rotateLeft(seedTab[kmerSeq.charAt(i)&cpOff], i);
        return hVal;
    }

    /**
     * ntHash basic function, i.e. ntBase, using rotate ops
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @return          hash value
     */
    public static long NT64(final CharSequence kmerSeq, final int k) {
        return getFhval(kmerSeq, k);
    }

    /**
     * ntHash with seeding option, using rotate ops
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param seed      seed value
     * @return          hash value
     */
    public static long NT64(final CharSequence kmerSeq, final int k, final int seed) {
        long hVal = getFhval(kmerSeq, k);
        hVal *= seed ^ k * multiSeed;
        hVal ^= hVal >>> multiShift;
        return hVal;
    }
    
    /**
     * ntHash for sliding k-mers, using rotate ops
     * @param fhVal     forward hash value
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @return          hash value
     */
    public static long NT64(final long fhVal, final char charOut, final char charIn, final int k) {
        return(Long.rotateLeft(fhVal, 1) ^ Long.rotateLeft(seedTab[charOut], k) ^ seedTab[charIn]);
    }

    /**
     * ntHash for sliding k-mers, using rotate ops
     * @param fhVal     forward hash value
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @return          hash value
     */
    public static long NT64B(final long fhVal, final char charOut, final char charIn, final int k) {
        return(Long.rotateRight(fhVal, 1) ^ Long.rotateRight(seedTab[charOut], 1) ^ Long.rotateLeft(seedTab[charIn], k-1));
    }
    
    /**
     * Canonical ntBase, using rotate ops
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @return          hash value
     */
    public static long NTC64(final CharSequence kmerSeq, final int k) {
        long fhVal = getFhval(kmerSeq, k);
        long rhVal = getRhval(kmerSeq, k);
        return (rhVal<fhVal)? rhVal : fhVal;
    }

    /**
     * Canonical ntHash for sliding k-mers, using rotate ops
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @param frhVals   array to store the forward and reverse hash values of kmer
     * @return          hash value
     */
    public static long NTC64(final char charOut, final char charIn, final int k, final long[] frhVals) {
        frhVals[0] = Long.rotateLeft(frhVals[0], 1) ^ Long.rotateLeft(seedTab[charOut], k) ^ seedTab[charIn]; // forward strand
        frhVals[1] = Long.rotateRight(frhVals[1], 1) ^ Long.rotateRight(seedTab[charOut&cpOff], 1) ^ Long.rotateLeft(seedTab[charIn&cpOff], k-1); // reverse strand
        return (frhVals[1]<frhVals[0])? frhVals[1] : frhVals[0]; // canonical
    }

    /**
     * Canonical ntHash for backward-sliding k-mers, using rotate ops
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @param frhVals   array to store the forward and reverse hash values of kmer
     * @return          hash value
     */
    public static long NTC64B(final char charOut, final char charIn, final int k, final long[] frhVals) {
        frhVals[0] = Long.rotateRight(frhVals[0], 1) ^ Long.rotateRight(seedTab[charOut], 1) ^ Long.rotateLeft(seedTab[charIn], k-1); // forward strand
        frhVals[1] = Long.rotateLeft(frhVals[1], 1) ^ Long.rotateLeft(seedTab[charOut&cpOff], k) ^ seedTab[charIn&cpOff]; // reverse strand
        return (frhVals[1]<frhVals[0])? frhVals[1] : frhVals[0]; // canonical
    }
    
    /**
     * Canonical ntBase with seeding option, using rotate ops
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param seed      seed value
     * @return          hash value
     */
    public static long NTC64(final CharSequence kmerSeq, final int k, final int seed) {
        long hVal = NTC64(kmerSeq,k);
        hVal *= seed ^ k * multiSeed;
        hVal ^= hVal >>> multiShift;
        return hVal;
    }

    /*
     * Using pre-computed seed matrix msTab instead of rotate opts
    */

    /**
     * ntBase
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @return          hash value
     */
    public static long NTP64(final CharSequence kmerSeq, final int k) {
        long hVal=0;
        for(int i=0; i<k; ++i)
            hVal ^= msTab[kmerSeq.charAt(i)][(k-1-i)%64];
        return hVal;
    }

    /**
     * ntBase (reverse complement)
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @return          hash value
     */    
    public static long NTP64RC(final CharSequence kmerSeq, final int k) {
        long rhVal=0;
        for(int i=0; i<k; ++i) {
            rhVal ^= msTab[kmerSeq.charAt(i)&cpOff][i%64];
        }
        return rhVal;
    }

    /**
     * ntBase with seeding option
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param seed      seed value
     * @return          hash value
     */
    public static long NTP64(final CharSequence kmerSeq, final int k, final int seed) {
        long hVal = NTP64(kmerSeq, k);
        hVal *= seed ^ k * multiSeed;
        hVal ^= hVal >>> multiShift;
        return hVal;
    }
    
    /**
     * ntHash for sliding k-mers
     * @param fhVal     forward hash value
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @return          hash value
     */
    public static long NTP64(final long fhVal, final char charOut, final char charIn, final int k) {
        return(Long.rotateLeft(fhVal, 1) ^ msTab[charOut][k%64] ^ msTab[charIn][0]);
    }

    /**
     * ntHash for backward-sliding k-mers
     * @param fhVal     forward hash value
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @return          hash value
     */
    public static long NTP64B(final long fhVal, final char charOut, final char charIn, final int k) {
        return(Long.rotateRight(fhVal, 1) ^ msTab[charOut][63] ^ msTab[charIn][(k-1)%64]);
    }
    
    /**
     * Canonical ntBase
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @return          hash value
     */
    public static long NTPC64(final CharSequence kmerSeq, final int k) {
        long fhVal=0, rhVal=0;
        for(int i=0; i<k; ++i) {
            fhVal ^= msTab[kmerSeq.charAt(i)][(k-1-i)%64];
            rhVal ^= msTab[kmerSeq.charAt(i)&cpOff][i%64];
        }
        return (rhVal<fhVal)? rhVal : fhVal;
    }
    
    /**
     * Canonical ntBase with seeding option
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param seed      seed value
     * @return          hash value
     */
    public static long NTPC64(final CharSequence kmerSeq, final int k, final int seed) {
        long hVal = NTPC64(kmerSeq, k);
        hVal *= seed ^ k * multiSeed;
        hVal ^= hVal >>> multiShift;
        return hVal;
    }
    
    /**
     * Canonical ntBase with forward/reverse hash values
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param frhVals   array to store the forward and reverse strand hash values of the kmer
     * @return          hash value
     */
    public static long NTPC64(final CharSequence kmerSeq, final int k, long[] frhVals) {
        frhVals[0]=0;
        frhVals[1]=0;
        for(int i=0; i<k; ++i) {
            frhVals[0] ^= msTab[kmerSeq.charAt(i)][(k-1-i)%64];
            frhVals[1] ^= msTab[kmerSeq.charAt(i)&cpOff][i%64];
        }
        return (frhVals[1]<frhVals[0])? frhVals[1] : frhVals[0];
    }
    

    /**
     * Canonical ntHash for sliding k-mers
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @param frhVals   array to store the forward and reverse strand hash values of the kmer
     * @return          hash value
     */
    public static long NTPC64(final char charOut, final char charIn, final int k, final long[] frhVals) {
        frhVals[0] = Long.rotateLeft(frhVals[0], 1) ^ msTab[charOut][k%64] ^ msTab[charIn][0]; // forward strand
        frhVals[1] = Long.rotateRight(frhVals[1], 1) ^ msTab[charOut&cpOff][63] ^ msTab[charIn&cpOff][(k-1)%64]; // reverse strand
        return (frhVals[1]<frhVals[0])? frhVals[1] : frhVals[0]; // canonical
    }

    /**
     * Canonical ntHash for backward-sliding k-mers
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @param frhVals   array to store the forward and reverse strand hash values of the kmer
     * @return          hash value
     */
    public static long NTPC64B(final char charOut, final char charIn, final int k, final long[] frhVals) {
        frhVals[0] = Long.rotateRight(frhVals[0], 1) ^ msTab[charOut][0] ^ msTab[charIn][k%64]; // forward strand
        frhVals[1] = Long.rotateLeft(frhVals[1], 1) ^ msTab[charOut&cpOff][(k-1)%64] ^ msTab[charIn&cpOff][63]; // reverse strand
        return (frhVals[1]<frhVals[0])? frhVals[1] : frhVals[0]; // canonical
    }
    
    /**
     * Multihash ntBase
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param hVal      array to store hash values
     */
    public static void NTM64(final CharSequence kmerSeq, final int k, final int m, final long[] hVal) {
        long bVal, tVal;
        bVal = NTP64(kmerSeq, k);
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }

    /**
     * Multihash ntBase (reverse complement)
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param hVal      array to store hash values
     */
    public static void NTM64RC(final CharSequence kmerSeq, final int k, final int m, final long[] hVal) {
        long bVal, tVal;
        bVal = NTP64RC(kmerSeq, k);
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }
    
    /**
     * Multihash ntHash for sliding k-mers
     * @param charOut   nucleotide to remove from the left
     * @param charIn    nucleotide to add to the right
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param hVal      array to store hash values
     */
    public static void NTM64(final char charOut, final char charIn, final int k, final int m, final long[] hVal) {
        long bVal, tVal;
        bVal = Long.rotateLeft(hVal[0], 1) ^ msTab[charOut][k%64] ^ msTab[charIn][0];
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }
    
    /**
     * Multihash ntHash for sliding k-mers
     * @param charOut   nucleotide to remove from the left
     * @param charIns    nucleotides to add to the right
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param hVal      array to store hash values
     * @param kMod64
     * @return hash values for each nucleotide added
     */
    public static long[][] NTM64(final char charOut, final char[] charIns, final int k, final int m, final long[] hVal, final int kMod64) {
        long bVal, tVal;
        
        final long tmpVal = Long.rotateLeft(hVal[0], 1) ^ msTab[charOut][kMod64];
        
        int numIns = charIns.length;
        final long[][] results = new long[numIns][m];
        
        for (int c=0; c<numIns; ++c) {
            bVal = tmpVal ^ msTab[charIns[c]][0];
            results[c][0] = bVal;
            for(int i=1; i<m; ++i) {
                tVal = bVal * (i ^ k * multiSeed);
                tVal ^= tVal >>> multiShift;
                results[c][i] = tVal;
            }
        }
        
        return results;
    }
    
    /**
     * Multihash ntHash for sliding reverse-complement k-mers
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param hVal      array to store hash values
     */
    public static void NTM64RC(final char charOut, final char charIn, final int k, final int m, final long[] hVal) {
        long bVal, tVal;
        bVal = Long.rotateRight(hVal[0], 1) ^ Long.rotateRight(seedTab[charOut&cpOff], 1) ^ Long.rotateLeft(seedTab[charIn&cpOff], k-1);
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }
    
    /**
     * Multihash ntHash for backward-sliding k-mers
     * @param charOut   nucleotide to remove from the right
     * @param charIn    nucleotide to add to the left
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param hVal      array to store hash values
     * @param kMinus1Mod64
     */
    public static void NTM64B(final char charOut, final char charIn, final int k, final int m, final long[] hVal, final int kMinus1Mod64) {
        long bVal, tVal;
        bVal = Long.rotateRight(hVal[0], 1) ^ msTab[charOut][63] ^ msTab[charIn][kMinus1Mod64];
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }
    
    /**
     * Multihash ntHash for backward-sliding k-mers
     * @param charOut   nucleotide to remove from the right
     * @param charIns    nucleotides to add to the left
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param hVal      array to store hash values
     * @param kMinus1Mod64
     * @return hash values for each nucleotide added
     */
    public static long[][] NTM64B(final char charOut, final char[] charIns, final int k, final int m, final long[] hVal, final int kMinus1Mod64) {
        long bVal, tVal;
        
        final long tmpVal = Long.rotateRight(hVal[0], 1) ^ msTab[charOut][63];
        
        int numIns = charIns.length;
        final long[][] results = new long[numIns][m];
        
        for (int c=0; c<numIns; ++c) {
            bVal = tmpVal ^ msTab[charIns[c]][kMinus1Mod64];
            results[c][0] = bVal;
            for(int i=1; i<m; ++i) {
                tVal = bVal * (i ^ k * multiSeed);
                tVal ^= tVal >>> multiShift;
                results[c][i] = tVal;
            }
        }
        
        return results;
    }
    
    /**
     * Canonical Multihash ntBase
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param hVal      array to store hash values
     */
    public static void NTMC64(final CharSequence kmerSeq, final int k, final int m, final long[] hVal) {
        long bVal, tVal;
        bVal = NTPC64(kmerSeq, k);
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }
 
    /**
     * Canonical Multihash ntHash
     * @param kmerSeq   kmer to be hashed
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param frhVals   array to store the forward and reverse strand hash values of the kmer
     * @param hVal      array to store hash values
     */
    public static void NTMC64(final CharSequence kmerSeq, final int k, final int m, final long[] frhVals, final long[] hVal) {
        long bVal, tVal;
        bVal = NTPC64(kmerSeq, k, frhVals);
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }

    /**
     * Canonical Multihash ntHash for sliding k-mers
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param frhVals   array to store the forward and reverse strand hash values of the kmer
     * @param hVal      array to store hash values
     */
    public static void NTMC64(final char charOut, final char charIn, final int k, final int m, final long[] frhVals, final long[] hVal) {
        long bVal, tVal;
        bVal = NTPC64(charOut, charIn, k, frhVals);
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }

    /**
     * Canonical Multihash ntHash for backward-sliding k-mers
     * @param charOut   nucleotide to remove
     * @param charIn    nucleotide to add
     * @param k         length of kmer
     * @param m         number of hash values to generate
     * @param frhVals   array to store the forward and reverse strand hash values of the kmer
     * @param hVal      array to store hash values
     */
    public static void NTMC64B(final char charOut, final char charIn, final int k, final int m, final long[] frhVals, final long[] hVal) {
        long bVal, tVal;
        bVal = NTPC64B(charOut, charIn, k, frhVals);
        hVal[0] = bVal;
        for(int i=1; i<m; ++i) {
            tVal = bVal * (i ^ k * multiSeed);
            tVal ^= tVal >>> multiShift;
            hVal[i] = tVal;
        }
    }
            
//    public static void main(String[] args) {
//        long h = Long.MIN_VALUE;
//        
//        System.out.println(Long.toBinaryString(h));
//                
//        System.out.println(Long.toBinaryString(Long.rotateLeft(h, 1)));
////        System.out.println(Long.toBinaryString(rol(h, 1)));
//        
//        System.out.println(Long.toBinaryString(Long.rotateRight(h, 1)));
////        System.out.println(Long.toBinaryString(ror(h, 1)));
//        
//        int k = 25;
//        int m = 2;
//
//        //            "123456789012345678901234567890"
//        String seq =  "AAAAAAAAAAAAAAAAAAAAAAAAA";
//        long hVal = NTP64(seq, k);
//        System.out.println(hVal);
//        
//        String seq0 = "CAAAAAAAAAAAAAAAAAAAAAAAA";
//        System.out.println(NTP64(seq0, k));
//
//
//        System.out.println(NTP64B(hVal, 'A', 'C', k));
//    }
}
