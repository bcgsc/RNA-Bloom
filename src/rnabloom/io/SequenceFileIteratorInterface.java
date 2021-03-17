/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.io;

import java.io.IOException;

/**
 *
 * @author Ka Ming Nip
 */
public interface SequenceFileIteratorInterface {
    public String next() throws IOException;
}
