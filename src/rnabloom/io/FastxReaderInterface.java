/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.IOException;

/**
 *
 * @author gengar
 */
public interface FastxReaderInterface {
    boolean hasNext();
    String next() throws FileFormatException;
    void close() throws IOException;
}
