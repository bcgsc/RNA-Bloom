/* 
 * Copyright (C) 2018-present BC Cancer Genome Sciences Centre
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
package rnabloom.io;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import static rnabloom.util.SeqBitsUtils.bitsToSeq;
import static rnabloom.util.SeqBitsUtils.seqLenToNumBytes;
import static rnabloom.util.SeqBitsUtils.fourBytesToInt;

/**
 *
 * @author Ka Ming Nip
 */
public class NucleotideBitsReader implements SequenceFileIteratorInterface {
    private final FileInputStream fin;
    private final byte[] seqLenBytes = new byte[4];
    
    public NucleotideBitsReader(String path) throws FileNotFoundException, IOException {
        fin = new FileInputStream(path);
    }
        
    @Override
    public synchronized String next() throws IOException {
        if (fin.read(seqLenBytes) > 0) {
            int seqLen = fourBytesToInt(seqLenBytes);
            byte[] seqBytes = new byte[seqLenToNumBytes(seqLen)];
            if (fin.read(seqBytes) > 0) {
                return bitsToSeq(seqBytes, seqLen);
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
            while ((seq=reader.next()) != null) {
                seqs2.add(seq);
            }
            reader.close();
            
            for (int i=0; i<seqs.length; ++i) {
                System.out.println(seqs[i].equals(seqs2.get(i)));
            }
            
        } catch (IOException ex) {
            System.out.println(ex.getMessage());
        }
        
    }
}
