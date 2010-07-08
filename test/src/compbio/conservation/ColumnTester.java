package compbio.conservation;

import org.testng.Assert;
import org.testng.annotations.Test;

public class ColumnTester {
	
	@Test
	public void OccuranceTester() {
		
		char a = 'D';
		char b = 'D';
		char c = 'D';
		char d = 'D';
		char e = 'D';
		char f = 'D';
		char g = 'D';
		char h = 'D';
		char i = 'D';
		char j = '-';
	
		
	AminoAcidColumn col = new AminoAcidColumn(a, b, c, d, e, f, g, h, i, j);
	
	int nr = col.getOccurance('D');
	Assert.assertEquals(nr, 9);
	int res = col.howManyResidueTypes();
	Assert.assertEquals(res, 1);
	Assert.assertEquals(col.allButOneGaps(), false);
	
	
	
		
	}

	@Test
	
	public void karlinScoreTester() {
		
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
	
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double kar = col.karlinScore();
	
	Assert.assertTrue(kar != -2.0);
	Assert.assertEquals(kar, 2.0, 0.1);
	
	}

	@Test
	public void armonScoreTester() {
		
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
	
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double arm = col.armonScore();
	
	Assert.assertEquals(arm, 2.0, 0.1);
		
	}
	
	@Test
	public void thompsonScoreTester() {
		
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
	
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double thom = col.thompsonScore();
	
	Assert.assertEquals(thom, 2.0, 0.1);
	
	}
		
	@Test
	public void lancetScoreTester() {
		
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
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double lan = col.notLancetScore();
	
	Assert.assertEquals(lan, 2.0, 0.1);
	
	}
	
	@Test
	public void mirnyTester() {
		
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
	
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double mirny = col.mirnyScore();
	
	Assert.assertEquals(mirny, 2.0, 0.1);
	
		
	}
		
	@Test
	public void WilliamsonTester() {
		
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
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double will = col.williamsonScore();
	
	Assert.assertEquals(will, 2.0, 0.1);
	
		
	}

	@Test
	public void TaylorTester() {
		
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
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);	
	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	int tay = col.SmallestTaylorSetGaps();
	Assert.assertEquals(tay, 3);
	
	}
	
	@Test
	public void KabatTester() {
		
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
	
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double kab = col.kabatScore();
	
	Assert.assertEquals(kab, 1.0, 0.1);
	
	}
	
	@Test
	public void schneiderTester() {
		
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
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);	
	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double sch = col.schneiderScore();
	
	Assert.assertEquals(sch, 1.0, 0.1);
	
	}
	
	@Test
	public void joresScoreTester() {
		
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
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);	
	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double jor = col.joresScore();
	
	Assert.assertEquals(jor, 2.0, 0.1);
	
	}
	
	
	@Test
	public void sanderScoreTester() {
		
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
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double san = col.sanderScore();
	
	Assert.assertEquals(san, 2.0, 0.1);
	
	}
	
	@Test
	public void valdarScoreTester() {
		
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
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double val = col.valdarScore();
	
	Assert.assertEquals(val, 2.0, 0.1);
	
	}
	
	@Test
	public void landgarfScoreTester() {
		
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
		
	AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
	
		
	Column col = new Column(a, b, c, d, e, f, g, h, i, j, matrix);
	
	double lan = col.landgrafScore();
	
	Assert.assertEquals(lan, 2.0, 0.1);
	
	}
	
	
}

