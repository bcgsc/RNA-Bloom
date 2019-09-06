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
public class PafRecord {
    public String qName = null, tName = null;
    public boolean reverseComplemented = false;
    public int qLen, tLen, qStart, qEnd, tStart, tEnd, numMatch = -1;

    public PafRecord() {
    }

    public PafRecord(String[] cols) {
        qName = cols[0];
        qLen = Integer.parseInt(cols[1]);
        qStart = Integer.parseInt(cols[2]);
        qEnd = Integer.parseInt(cols[3]);
        reverseComplemented = cols[4].equals("-");
        tName = cols[5];
        tLen = Integer.parseInt(cols[6]);
        tStart = Integer.parseInt(cols[7]);
        tEnd = Integer.parseInt(cols[8]);
        numMatch = Integer.parseInt(cols[9]);
    }

    public void update(String[] cols) {
        qName = cols[0];
        qLen = Integer.parseInt(cols[1]);
        qStart = Integer.parseInt(cols[2]);
        qEnd = Integer.parseInt(cols[3]);
        reverseComplemented = cols[4].equals("-");
        tName = cols[5];
        tLen = Integer.parseInt(cols[6]);
        tStart = Integer.parseInt(cols[7]);
        tEnd = Integer.parseInt(cols[8]);
        numMatch = Integer.parseInt(cols[9]);
    }
}
