/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package matcher;

/**
 *
 * @author jiadonglin
 */
public class TestStringMatcher {

    public static void main(String[] args) {
        StringMatcher matcher = new StringMatcher();

//        String ref = "TTGTTTGTCGAGAGCAG";
        String ref = "ATGTTATTTTTACAATGCTAAAAATTAGAGTAGAAAATTAAAAAGTAAAGTCTAGATTGTTCAATCGGGCTTGTTTTAAATTTGTTTTTTATATTTTTCACTCTAGAGAAGAGCACCTAGCATCTTAACTCCTCATGCTGAAGAGGAAATAGTTCTGAGCAATTATTTATTTAAAATATTAAATAATATTTAAAATATTAAATAATGTTTAAAATATTTAAAATATTAAATAATGTTTAAAATATTTAAAATATTAAATAATGTTTAAAATATTTAAAATATTAAATAATGTTTAAAATATTTAAAATATTAATGTTTAAAATATTTAAAATATTAAATAATGTTTAAAATATTTAAAATATTAAATAATATTTAAAATATTTAAAACATTAGTGGTGGT";
        String read = "TTGTTTTTTATATTTTTCACTCTAGAGAATAGCACCTAGCATCTTAACTCCTCATGCTGAAGAGGAAATAGTTCTGAGCAATTATTTATTTAAAATATTAAATAATATTTAAAATATTTAAAATATTAAAT";


//        String read = "TTGTTAGAGAAG";
//        System.out.println(ref.substring(100, 178));

        System.out.println("Ref len: " + ref.length() + " Read len: " + read.length());
        matcher.alignStrToRef(read, ref);


    }
}
