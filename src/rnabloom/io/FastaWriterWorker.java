/* 
 * Copyright (C) 2021-present BC Cancer Genome Sciences Centre
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

import java.util.Queue;

/**
 *
 * @author Ka Ming Nip
 */
public class FastaWriterWorker implements Runnable {
    private final Queue<String> inputQueue;
    private final FastaWriter writer;
    private boolean terminateWhenInputExhausts = false;
    private boolean successful = false;
    private Exception exception = null;
    private int numSeq = 0;
    private final String seqPrefix;

    public FastaWriterWorker(Queue<String> inputQueue, 
            FastaWriter writer, String seqPrefix) {
        this.inputQueue = inputQueue;
        this.writer = writer;
        this.seqPrefix = seqPrefix;
    }

    @Override
    public void run() {
        try {
            while(true) {
                String seq = inputQueue.poll();

                if (seq == null) {
                    if (terminateWhenInputExhausts) {
                        break;
                    }
                    continue;
                }
                
                String header = seqPrefix + ++numSeq + " l=" + Integer.toString(seq.length());

                writer.write(header, seq);
            }

            successful = true;

        } catch (Exception ex) {
            System.out.println(ex.getMessage());
            exception = ex;
            successful = false;
        }
    }

    public void terminateWhenInputExhausts() {
        terminateWhenInputExhausts = true;
    }

    public boolean isSucessful() {
        return successful;
    }

    public Exception getExceptionCaught() {
        return exception;
    }

    public int getNumSequencesWritten() {
        return numSeq;
    }
}
