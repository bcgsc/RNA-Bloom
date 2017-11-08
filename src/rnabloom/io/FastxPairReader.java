/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.IOException;

/**
 *
 * @author kmnip
 */
public interface FastxPairReader {
    
    public boolean hasNext();
    
    public PairedReadSegments next() throws FileFormatException;
    
    public void close() throws IOException;
}
