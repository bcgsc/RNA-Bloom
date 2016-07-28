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
abstract class AbstractBloomFilter{
    abstract void add(String key);
    abstract boolean lookup(String key);
    abstract void save(String path);
    abstract void restore(String path);
    abstract float getFPR();
}
