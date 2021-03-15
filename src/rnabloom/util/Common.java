/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package rnabloom.util;

import java.math.BigDecimal;
import java.math.MathContext;

/**
 *
 * @author Ka Ming Nip
 */
public class Common {
    public static float convertToRoundedPercent(float f) {
        return roundToSigFigs(f * 100, 3);
    }
    
    public static float roundToSigFigs(float f, int sigFigs) {
        BigDecimal bd = new BigDecimal(f);
        bd = bd.round(new MathContext(sigFigs));
        return bd.floatValue();
    }
}
