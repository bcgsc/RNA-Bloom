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
public class FastxFilePair {
    public final String leftPath;
    public final String rightPath;
    public final boolean leftRevComp;
    public final boolean rightRevComp;

    public FastxFilePair(String leftPath, String rightPath, boolean leftRevComp, boolean rightRevComp) {
        this.leftPath = leftPath;
        this.rightPath = rightPath;
        this.leftRevComp = leftRevComp;
        this.rightRevComp = rightRevComp;
    }
}
