/* 
 * Copyright (C) 2018-present BC Cancer Genome Sciences Centre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
