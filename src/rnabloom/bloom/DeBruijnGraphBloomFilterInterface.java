/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

/**
 *
 * @author kmnip
 */
interface DeBruijnGraphBloomFilterInterface extends BloomFilterInterface, CountingBloomFilterInterface{
        
    public void add(String kmer);
    public boolean lookup(String kmer);
    public void increment(String kmer);
    public float getCount(String kmer);
    public float getFPR();
    public String[] getPredecessors(String kmer);
    public String[] getSuccessors(String kmer);
    
}
