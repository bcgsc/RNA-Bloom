/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

/**
 *
 * @author kmnip
 */
public class FastqRecord {
    public String name;
    public String seq;
    public String qual;    
    
    public FastqRecord() {
        
    }
    
    public FastqRecord(String name, String seq, String qual) {
        this.name = name;
        this.seq = seq;
        this.qual = qual;
    }
}
