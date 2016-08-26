/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.function.Supplier;
import java.util.regex.Pattern;
import rnabloom.io.FastqPairReader.ReadPair;
import static rnabloom.util.SeqUtils.filterFastq;
import static rnabloom.util.SeqUtils.reverseComplement;

/**
 *
 * @author kmnip
 */
public class FastqPairReader implements Iterator<ReadPair> {
    private final FastqReader leftReader;
    private final FastqReader rightReader;
    private final Pattern qualPattern;
    private final Supplier<ReadPair> nextFunction;
                
    public FastqPairReader(FastqReader leftReader, FastqReader rightReader, Pattern qualPattern, boolean revCompLeft, boolean revCompRight) {
        this.qualPattern = qualPattern;
        this.leftReader = leftReader;
        this.rightReader = rightReader;
        if (revCompLeft) { 
            if (revCompRight) {
                nextFunction = this::nextRR;
            }
            else {
                nextFunction = this::nextRF;
            }
        }
        else {
            if (revCompRight) {
                nextFunction = this::nextFR;
            }
            else {
                nextFunction = this::nextFF;
            }                
        }
    }
    
    public static final class ReadPair {
        public String left;
        public String right;
    }

    @Override
    public boolean hasNext() {
        return this.leftReader.hasNext() && this.rightReader.hasNext();
    }

    @Override
    public ReadPair next() {
        return nextFunction.get();
    }

    private ReadPair nextFF() {
        ReadPair p = new ReadPair();

        FastqRecord frLeft = leftReader.next();
        FastqRecord frRight = rightReader.next();

        if (!frLeft.name.equals(frRight.name)) {
            throw new NoSuchElementException("Inconsistent record names: \"" + frLeft.name + "\" and \"" + frRight.name + "\"");
        }

        ArrayList<String> lefts = filterFastq(frLeft, qualPattern);
        if (lefts.isEmpty()) {
            p.left = "";
        }
        else if (lefts.size() == 1) {
            p.left = lefts.get(0);
        }
        else {
            p.left = lefts.get(0);
            int longest = p.left.length();
            for (String l : lefts) {
                int len = l.length();
                if (len > longest) {
                    p.left = l;
                    longest = len;
                }
            }
        }

        ArrayList<String> rights = filterFastq(frRight, qualPattern);
        if (rights.isEmpty()) {
            p.right = "";
        }
        else if (lefts.size() == 1) {
            p.right = rights.get(0);
        }
        else {
            p.right = rights.get(0);
            int longest = p.right.length();
            for (String r : rights) {
                int len = r.length();
                if (len > longest) {
                    p.left = r;
                    longest = len;
                }
            }
        }

        return p;
    }

    private ReadPair nextFR() {
        ReadPair p = nextFF();
        p.right = reverseComplement(p.right);
        return p;
    }

    private ReadPair nextRF() {
        ReadPair p = nextFF();
        p.left = reverseComplement(p.left);
        return p;
    }

    private ReadPair nextRR() {
        ReadPair p = nextFF();
        p.left = reverseComplement(p.left);
        p.right = reverseComplement(p.right);
        return p;
    } 
}
