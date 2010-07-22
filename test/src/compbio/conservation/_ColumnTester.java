package compbio.conservation;

import org.testng.Assert;
import org.testng.annotations.Test;

public class _ColumnTester {
	
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
	
	_Column col = new _Column(column, matrix);
	

	@Test
	
	public void karlinScoreTester() {
		
	double kar = col.karlinScore();
	
	Assert.assertTrue(kar != -2.0);
	Assert.assertEquals(kar, 0.269, 0.001);
	
	}

	@Test
	public void armonScoreTester() {
	
	double arm = col.armonScore();
	
	Assert.assertEquals(arm, 26.760, 0.001);
		
	}
	
	@Test
	public void thompsonScoreTester() {
		
	double thom = col.thompsonScore();
	
	Assert.assertEquals(thom, 8.784, 0.001);
	
	}
		
	@Test
	public void lancetScoreTester() {
	
	double lan = col.notLancetScore();
	
	Assert.assertEquals(lan, -0.490, 0.001);
	
	}
	
	@Test
	public void mirnyTester() {
		
	double mirny = col.mirnyScore();
	
	Assert.assertEquals(mirny, -0.409, 0.001);
	
		
	}
		
	@Test
	public void WilliamsonTester() {

	double will = col.williamsonScore();
	
	Assert.assertEquals(will, -2.072, 0.001);
	
		
	}

	@Test
	public void TaylorTester() {

	int tay = col.SmallestTaylorSetGaps();
	Assert.assertEquals(tay, 21);
	
	}
	
	@Test
	public void KabatTester() {
		
	double kab = col.kabatScore();
	
	Assert.assertEquals(kab, 6.0, 0.001);
	
	}
	
	@Test
	public void schneiderTester() {

	double sch = col.schneiderScore();
	
	Assert.assertEquals(sch, 0.313, 0.001);
	
	}
	
	@Test
	public void joresScoreTester() {

	double jor = col.joresScore();
	
	Assert.assertEquals(jor, 15.0, 0.001);
	
	}
	
	
	@Test
	public void sanderScoreTester() {

	double san = col.sanderScore();
	
	Assert.assertEquals(san, -900.000, 0.001);
	
	}
	
	@Test
	public void valdarScoreTester() {

	double val = col.valdarScore();
	
	Assert.assertEquals(val, 7.896, 0.001);
	
	}
	
	@Test
	public void landgarfScoreTester() {

	double lan = col.landgrafScore();
	
	//Assert.assertEquals(lan, 1177.130, 0.001);
	
	}
	
	
}

