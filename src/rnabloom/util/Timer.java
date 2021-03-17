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

/**
 *
 * @author Ka Ming Nip
 */
public class Timer {
    private final long globalStartTime;
    private long startTime;

    public Timer() {
        globalStartTime = System.currentTimeMillis();
        startTime = globalStartTime;
    }

    public void start() {
        startTime = System.currentTimeMillis();
    }

    public long elapsedMillis() {
        return System.currentTimeMillis() - startTime;
    }

    public long totalElapsedMillis() {
        return System.currentTimeMillis() - globalStartTime;
    }

    public String elapsedDHMS() {
        return dhmsFormat(elapsedMillis());
    }

    public String totalElapsedDHMS() {
        return dhmsFormat(totalElapsedMillis());
    }
    
    public static String dhmsFormat(long millis) {
        long seconds = millis / 1000;

        long days = seconds / 86400;
        seconds = seconds % 86400;
        
        long hours = seconds / 3600;
        seconds = seconds % 3600;

        long minutes = seconds / 60;
        seconds = seconds % 60;

        StringBuilder sb = new StringBuilder();

        if (days > 0) {
            sb.append(days);
            sb.append("d ");
        }
        
        if (days > 0 || hours > 0) {
            sb.append(hours);
            sb.append("h ");
        }

        if (days > 0 || hours > 0 || minutes > 0) {
            sb.append(minutes);
            sb.append("m ");
        }

        if (days == 0 && hours == 0 && minutes == 0) {
            sb.append(millis/1000f);
            sb.append("s");
        }
        else {
            sb.append(seconds);
            sb.append("s");
        }

        return sb.toString();
    }
}
