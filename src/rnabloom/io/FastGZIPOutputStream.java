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
package rnabloom.io;

import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.Deflater;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author Ka Ming Nip
 */
public class FastGZIPOutputStream extends GZIPOutputStream {

    public FastGZIPOutputStream(OutputStream out) throws IOException {
        super(out);
        def.setLevel(Deflater.BEST_SPEED);
    }
    
    public FastGZIPOutputStream(OutputStream out, int size) throws IOException {
        super(out, size);
        def.setLevel(Deflater.BEST_SPEED);
    }

    public FastGZIPOutputStream(OutputStream out, boolean syncFlush) throws IOException {
        super(out, syncFlush);
        def.setLevel(Deflater.BEST_SPEED);
    }

    public FastGZIPOutputStream(OutputStream out, int size, boolean syncFlush) throws IOException {
        super(out, size, syncFlush);
        def.setLevel(Deflater.BEST_SPEED);
    }
    
}
