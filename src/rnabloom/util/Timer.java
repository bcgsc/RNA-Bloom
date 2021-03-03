/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
