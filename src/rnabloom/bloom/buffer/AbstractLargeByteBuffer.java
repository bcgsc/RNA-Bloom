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
package rnabloom.bloom.buffer;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 *
 * @author Ka Ming Nip
 */
public abstract class AbstractLargeByteBuffer {
    public abstract void set(long index, byte value);
    public abstract byte compareAndSwap(long index, byte expected, byte updated);
    public abstract byte get(long index);
    public abstract long size();
    public abstract long popCount();
    public abstract void empty();
    public abstract void destroy();
    public abstract void write(FileOutputStream out) throws IOException;
    public abstract void read(FileInputStream in) throws IOException;
}
