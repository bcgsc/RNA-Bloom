/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author kmnip
 */
public class NTCardHistogram {
    public long numKmers = 0;
    public int covThreshold = 0;
    
    public NTCardHistogram(String path) throws FileNotFoundException, IOException {
        boolean f0Found = false;
        long lastCount = -1;
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line;
        while ((line = br.readLine()) != null) {
            if (line.length() > 0) {
                String[] cols = line.split("\t");
                if (cols[0].equals("F0")) {
                    numKmers = Long.parseLong(cols[1]); 
                    f0Found = true;
                }
                else if (f0Found) {
                    long count = Long.parseLong(cols[1]);
                    if (lastCount <= 0) {
                        lastCount = count;
                    }
                    else if (3 * count > lastCount) {
                        covThreshold = Integer.parseInt(cols[0]) - 1;
                        break;
                    }
                    else {
                        lastCount = count;
                    }
                }
            }
        }
        br.close();
    }
    
}
