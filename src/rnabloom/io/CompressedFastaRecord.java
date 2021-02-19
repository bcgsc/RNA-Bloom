/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import rnabloom.util.BitSequence;

/**
 *
 * @author gengar
 */
public class CompressedFastaRecord {
    public String name;
    public BitSequence seqbits;
    
    public CompressedFastaRecord(String name, String seq) {
        this.name = name;
        this.seqbits = new BitSequence(seq); 
    }
}
