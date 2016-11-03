
import java.nio.ByteBuffer;
import java.nio.LongBuffer;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author gengar
 */
public class Tests {
        
    public static void main(String[] args) {
        long l = Long.MAX_VALUE;
        l &= ~(1 << 7);
        String lstr = Long.toBinaryString(l);
        //System.out.println(lstr);
        System.out.println(("0000000000000000000000000000000000000000000000000000000000000000" + lstr).substring(lstr.length()));
        
        byte b = Byte.MIN_VALUE;
        b |= (1 << 5);
        b |= (1 << 2);
        String bstr = Integer.toBinaryString(b);
        //System.out.println(bstr);
        System.out.println(("00000000" + bstr).substring(bstr.length()));
        
        long l2 = Long.MIN_VALUE;
        String lstr2 = Long.toBinaryString(l2);
        System.out.println(("0000000000000000000000000000000000000000000000000000000000000000" + lstr2).substring(lstr2.length()));
        
        /*
        ByteBuffer bb = ByteBuffer.allocateDirect(Long.BYTES);
        bb.asLongBuffer().put(l);
        bb.put(0, b);
        long l2 = bb.getLong(0);
        String lstr2 = Long.toBinaryString(l2);
        System.out.println(("0000000000000000000000000000000000000000000000000000000000000000" + lstr2).substring(lstr2.length()));
        */
    }
    
}
