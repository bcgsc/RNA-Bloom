/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.bloom;

/**
 *
 * @author kmnip
 */
interface BloomFilterInterface{
    void add(String key);
    boolean lookup(String key);
    float getFPR();
}
