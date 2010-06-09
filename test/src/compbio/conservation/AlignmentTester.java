package compbio.conservation;

import org.testng.Assert;
import org.testng.annotations.Test;
import static org.testng.Assert.assertTrue;

public class AlignmentTester {
	

		
	
	
	@Test
	public void testKabat() {
		
		char p1 = 'D'; char p2 = 'D'; char p3 = 'D'; char p4 = 'D'; char p5 = 'D'; char p6 = 'D'; char p7 = 'I'; char p8 = 'P'; char p9 = 'D'; char p10 = 'L';
        char p11 = 'D'; char p12 = 'D'; char p13 = 'D'; char p14 = 'D'; char p15 = 'D'; char p16 = 'D'; char p17 = 'I'; char p18 = 'P';char p19 = 'V'; char p20 = 'L';
        char p21 = 'D'; char p22 = 'D'; char p23 = 'D'; char p24 = 'D'; char p25 = 'D'; char p26 = 'D'; char p27 = 'I'; char p28 = 'P';char p29 = 'Y'; char p30 = 'L';
        char p31 = 'D'; char p32 = 'D'; char p33 = 'D'; char p34 = 'D'; char p35 = 'D'; char p36 = 'D'; char p37 = 'I'; char p38 = 'P';char p39 = 'A'; char p40 = 'L';
        char p41 = 'D'; char p42 = 'D'; char p43 = 'D'; char p44 = 'D'; char p45 = 'D'; char p46 = 'D'; char p47 = 'L'; char p48 = 'W';char p49 = 'T'; char p50 = '-';
        char p51 = 'D'; char p52 = 'D'; char p53 = 'E'; char p54 = 'D'; char p55 = 'E'; char p56 = 'E'; char p57 = 'L'; char p58 = 'W';char p59 = 'K'; char p60 = '-';
        char p61 = 'D'; char p62 = 'D'; char p63 = 'E'; char p64 = 'D'; char p65 = 'E'; char p66 = 'E'; char p67 = 'L'; char p68 = 'W';char p69 = 'P'; char p70 = '-';
        char p71 = 'D'; char p72 = 'D'; char p73 = 'E'; char p74 = 'D'; char p75 = 'E'; char p76 = 'E'; char p77 = 'L'; char p78 = 'W';char p79 = 'C'; char p80 = '-';
        char p81 = 'D'; char p82 = 'D'; char p83 = 'E'; char p84 = 'D'; char p85 = 'E'; char p86 = 'F'; char p87 = 'V'; char p88 = 'S';char p89 = 'R'; char p90= '-';
        char p91 = 'D'; char p92 = 'E'; char p93 = 'E'; char p94 = 'F'; char p95 = 'F'; char p96 = 'F'; char p97 = 'V'; char p98 = 'S';char p99 = 'H'; char p100  = '-';	
		
	final AminoAcidMatrix m = new AminoAcidMatrix( p1,  p2,  p3, p4, p5,  p6, p7 ,  p8 , p9 , p10, p11,  p12, p13, p14, p15, p16, p17 , p18 , p19 , p20,p21, p22, p23, p24, p25, p26, p27 ,p28 ,p29 ,p30, p31, p32, p33, p34, p35, p36, p37, p38 , p39 , p40, p41, p42, p43, p44, p45, p46, p47 , p48 , p49 , p50, p51, p52, p53, p54, p55, p56, p57 , p58 , p59 , p60, p61, p62, p63, p64, p65, p66, p67 , p68 , p69 , p70, p71, p72, p73, p74, p75, p76, p77 , p78 , p79 , p80, p81, p82, p83, p84, p85, p86, p87 , p88 , p89 , p90, p91, p92, p93, p94, p95, p96, p97 , p98 , p99 , p100);
		
		Alignment a = new Alignment(m);
		double[] result = a.kabat();
		assertTrue(result[1] != 0);

	}
	
	
	@Test
	public void testJores() {
		
		char p1 = 'D'; char p2 = 'D'; char p3 = 'D'; char p4 = 'D'; char p5 = 'D'; char p6 = 'D'; char p7 = 'I'; char p8 = 'P'; char p9 = 'D'; char p10 = 'L';
        char p11 = 'D'; char p12 = 'D'; char p13 = 'D'; char p14 = 'D'; char p15 = 'D'; char p16 = 'D'; char p17 = 'I'; char p18 = 'P';char p19 = 'V'; char p20 = 'L';
        char p21 = 'D'; char p22 = 'D'; char p23 = 'D'; char p24 = 'D'; char p25 = 'D'; char p26 = 'D'; char p27 = 'I'; char p28 = 'P';char p29 = 'Y'; char p30 = 'L';
        char p31 = 'D'; char p32 = 'D'; char p33 = 'D'; char p34 = 'D'; char p35 = 'D'; char p36 = 'D'; char p37 = 'I'; char p38 = 'P';char p39 = 'A'; char p40 = 'L';
        char p41 = 'D'; char p42 = 'D'; char p43 = 'D'; char p44 = 'D'; char p45 = 'D'; char p46 = 'D'; char p47 = 'L'; char p48 = 'W';char p49 = 'T'; char p50 = '-';
        char p51 = 'D'; char p52 = 'D'; char p53 = 'E'; char p54 = 'D'; char p55 = 'E'; char p56 = 'E'; char p57 = 'L'; char p58 = 'W';char p59 = 'K'; char p60 = '-';
        char p61 = 'D'; char p62 = 'D'; char p63 = 'E'; char p64 = 'D'; char p65 = 'E'; char p66 = 'E'; char p67 = 'L'; char p68 = 'W';char p69 = 'P'; char p70 = '-';
        char p71 = 'D'; char p72 = 'D'; char p73 = 'E'; char p74 = 'D'; char p75 = 'E'; char p76 = 'E'; char p77 = 'L'; char p78 = 'W';char p79 = 'C'; char p80 = '-';
        char p81 = 'D'; char p82 = 'D'; char p83 = 'E'; char p84 = 'D'; char p85 = 'E'; char p86 = 'F'; char p87 = 'V'; char p88 = 'S';char p89 = 'R'; char p90= '-';
        char p91 = 'D'; char p92 = 'E'; char p93 = 'E'; char p94 = 'F'; char p95 = 'F'; char p96 = 'F'; char p97 = 'V'; char p98 = 'S';char p99 = 'H'; char p100  = '-';	
		
	final AminoAcidMatrix m = new AminoAcidMatrix( p1,  p2,  p3, p4, p5,  p6, p7 ,  p8 , p9 , p10, p11,  p12, p13, p14, p15, p16, p17 , p18 , p19 , p20,p21, p22, p23, p24, p25, p26, p27 ,p28 ,p29 ,p30, p31, p32, p33, p34, p35, p36, p37, p38 , p39 , p40, p41, p42, p43, p44, p45, p46, p47 , p48 , p49 , p50, p51, p52, p53, p54, p55, p56, p57 , p58 , p59 , p60, p61, p62, p63, p64, p65, p66, p67 , p68 , p69 , p70, p71, p72, p73, p74, p75, p76, p77 , p78 , p79 , p80, p81, p82, p83, p84, p85, p86, p87 , p88 , p89 , p90, p91, p92, p93, p94, p95, p96, p97 , p98 , p99 , p100);
		
		Alignment a = new Alignment(m);
		double[] result = a.jores();
		Assert.assertEquals(result[0], 1.0);
		Assert.assertEquals(result[4], 2.0);
		
	}

	@Test
	public void testJores2() { 
		
		char a = 'D';
		char b = 'D';
		char c = 'D';
		char d = 'D';
		char e = 'D';
		char f = 'D';
		char g = 'D';
		char h = 'D';
		char i = 'D';
		char j = 'E';
	
		
		AminoAcidColumn[] z = new AminoAcidColumn[1];
		
		z[0] = new AminoAcidColumn(a, b, c, d, e, f, g, h, i, j);
		
		double[] result = Alignment.jores2(z);
		
		Assert.assertEquals(result[0], 2);
		
		
		
		
	}
	
}