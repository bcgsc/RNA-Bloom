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
public class FastaRecord {
    public String name;
    public String seq; 
    
    public FastaRecord() {
        
    }
    
    public FastaRecord(String name, String seq) {
        this.name = name;
        this.seq = seq;
    }
}
