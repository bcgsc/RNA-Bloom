/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import static rnabloom.util.NucleotideBitsUtils.byteArrayToSeq;

/**
 *
 * @author kmnip
 */
public class NucleotideBitsReader {
    private FileInputStream fin;
    private byte[] seqLenBytes = new byte[4];
    
    public NucleotideBitsReader(String path) throws FileNotFoundException, IOException {
        fin = new FileInputStream(path);
    }
        
    public String read() throws IOException {
        if (fin.read(seqLenBytes) > 0) {
            int seqLen = ByteBuffer.wrap(Arrays.copyOfRange(seqLenBytes, 0, 4)).getInt();
            byte[] seqBytes = new byte[seqLen % 4 > 0 ? seqLen/4+1 : seqLen/4];
            if (fin.read(seqBytes) > 0) {
                return byteArrayToSeq(seqLen, seqBytes);
            }
        }
        
        return null;
    }
    
    public void close() throws IOException {
        fin.close();
    }
    
    public static void main(String[] args) {
        String[] seqs = new String[]{"CACGAGACCTCTCTACATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAGGCAGCT",
                                    "CCCAGATGGGTCTTGTCCCAGGTGCGGCTACAGCAGTGGGG",
                                    "CGCAGGACTGTTGAAGCCTTCGGAGACCCTGTCTCTTATACAAATCTCCGAGCCCACGAGAGCTCTCAGGATA",
                                    "TCGTATGGATGAGACAGCTTGAACACACAA"};
        
        String path = "~/test.2bit";
        try {
            NucleotideBitsWriter writer = new NucleotideBitsWriter(path, false);
            for (String s : seqs) {
                writer.write(s);
            }
            writer.close();
            
            NucleotideBitsReader reader = new NucleotideBitsReader(path);
            java.util.ArrayList<String> seqs2 = new java.util.ArrayList<>();
            String seq;
            while ((seq=reader.read()) != null) {
                seqs2.add(seq);
            }
            reader.close();
            
            for (int i=0; i<seqs.length; ++i) {
                System.out.println(seqs[i].equals(seqs2.get(i)));
            }
            
        } catch (IOException ex) {
            Logger.getLogger(NucleotideBitsReader.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
}
