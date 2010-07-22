package compbio.conservation;

import org.testng.Assert;
import org.testng.annotations.Test;

public class ColumnScoresTester {

	char a = 'D';
	char b = 'D';
	char c = 'D';
	char d = 'D';
	char e = 'D';
	char f = 'E';
	char g = 'E';
	char h = 'E';
	char i = 'F';
	char j = '-';
	
	char[] column = { a, b, c, d, e, f, g, h, i, j };
	
	AminoAcidMatrix matrix = new AminoAcidMatrix(column);
	

	@Test
	
	public void karlinScoreTester() {
		
	double kar = ColumnScores.karlinScore(matrix, 0);
	
	System.out.println(kar);
	
	//Assert.assertTrue(kar != -2.0);
	//Assert.assertEquals(kar, 0.269, 0.001);
	
	}

	@Test
	public void armonScoreTester() {
	
	double arm = ColumnScores.armonScore(matrix,0);
	
	System.out.println(arm);
	
	//Assert.assertEquals(arm, 26.760, 0.001);
		
	}
	
	@Test
	public void thompsonScoreTester() {
		
	double thom = ColumnScores.thompsonScore(matrix,0);
	
	System.out.println(thom);
	
	//Assert.assertEquals(thom, 8.784, 0.001);
	
	}
		
	@Test
	public void lancetScoreTester() {
	
	double lan = ColumnScores.notLancetScore(matrix, 0);
	
	System.out.println(lan);
	
	//Assert.assertEquals(lan, -0.490, 0.001);
	
	}
	
	@Test
	public void mirnyTester() {
		
	double mirny = ColumnScores.mirnyScore(matrix, 0);
	
	System.out.println(mirny);
	
	//Assert.assertEquals(mirny, -0.409, 0.001);
	
		
	}
		
	@Test
	public void WilliamsonTester() {

	double will = ColumnScores.williamsonScore(matrix, 0);
	
	System.out.println(will);
	
	//Assert.assertEquals(will, -2.072, 0.001);
	
		
	}

	@Test
	public void TaylorTester() {

	double tay = ColumnScores.taylorScoreGaps(matrix, 0);
	
	System.out.println(tay);
	
	//Assert.assertEquals(tay, 21);
	
	}
	
	@Test
	public void KabatTester() {
		
	double kab = ColumnScores.kabatScore(matrix, 0);
	
	System.out.println(kab);
	
	//Assert.assertEquals(kab, 6.0, 0.001);
	
	}
	
	@Test
	public void schneiderTester() {

	double sch = ColumnScores.schneiderScore(matrix, 0);
	
	System.out.println(sch);
	
	//Assert.assertEquals(sch, 0.313, 0.001);
	
	}
	
	@Test
	public void joresScoreTester() {

	double jor = ColumnScores.joresScore(matrix, 0);
	
	System.out.println(jor);
	
	//Assert.assertEquals(jor, 15.0, 0.001);
	
	}
	
	
	@Test
	public void sanderScoreTester() {

	double san = ColumnScores.sanderScore(matrix, 0);
	
	System.out.println(san);
	
	//Assert.assertEquals(san, -900.000, 0.001);
	
	}
	
	@Test
	public void valdarScoreTester() {

	double val = ColumnScores.valdarScore(matrix, 0);
	
	System.out.println(val);
	
	//Assert.assertEquals(val, 7.896, 0.001);
	
	}
	
	@Test
	public void landgarfScoreTester() {

	double lan = ColumnScores.landgrafScore(matrix, 0);
	
	System.out.println("Landgraf");
	
	System.out.println(lan);
	
	//Assert.assertEquals(lan, 1177.130, 0.001);
	
	}
	
	
}
