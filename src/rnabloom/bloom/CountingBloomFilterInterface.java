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
interface CountingBloomFilterInterface {
    void increment(String key);
    float getCount(String key);
    float getFPR();
}
