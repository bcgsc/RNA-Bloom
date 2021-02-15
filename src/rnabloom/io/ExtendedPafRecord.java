/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author kmnip
 */
public class ExtendedPafRecord extends PafRecord {
    public String cigar = null; // CIGAR string, eg. "129M" in "cg:Z:129M"
    public int nm = -1; // Total number of mismatches and gaps in the alignment, eg. 2 in "NM:i:2"
    public final static Pattern CIGAR_PATTERN = Pattern.compile("cg:Z:(\\S+)");
    public final static Pattern NM_PATTERN = Pattern.compile("NM:i:(\\d+)");
    private final Matcher cigarMatcher = CIGAR_PATTERN.matcher("");
    private final Matcher nmMatcher = NM_PATTERN.matcher("");
    
    @Override
    public void update(String[] cols) {
        super.update(cols);
        
        int numCols = cols.length;
        for (int i=12; i<numCols; ++i) {
            String item = cols[i];
            
            if (nm < 0) {
                nmMatcher.reset(item);
                if (nmMatcher.matches()) {
                    nm = Integer.parseInt(nmMatcher.group(1));

                    if (cigar != null) {
                        break;
                    }

                    continue;
                }
            }
            
            if (cigar == null) {
                cigarMatcher.reset(item);
                if (cigarMatcher.matches()) {
                    cigar = cigarMatcher.group(1);

                    if (nm >= 0) {
                        break;
                    }
                }
            }
        }
    }
}
