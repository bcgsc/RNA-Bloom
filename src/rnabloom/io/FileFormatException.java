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
public class FileFormatException extends IOException {
    
    public FileFormatException() {
        super();
    }
    
    public FileFormatException(String message) {
        super(message);
    }
    
    public FileFormatException(Exception e) {
        super(e);
    }
}
