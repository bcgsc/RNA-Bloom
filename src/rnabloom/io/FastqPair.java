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
public class FastqPair {
    public final String leftFastq;
    public final String rightFastq;
    public final boolean leftRevComp;
    public final boolean rightRevComp;

    public FastqPair(String leftFastq, String rightFastq, boolean leftRevComp, boolean rightRevComp) {
        this.leftFastq = leftFastq;
        this.rightFastq = rightFastq;
        this.leftRevComp = leftRevComp;
        this.rightRevComp = rightRevComp;
    }
}
