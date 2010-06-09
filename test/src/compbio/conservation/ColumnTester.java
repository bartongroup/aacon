package compbio.conservation;

import org.testng.Assert;
import org.testng.annotations.Test;

public class ColumnTester {
	
	@Test
	public void OccuranceTester() {
		
		char a = 'D';
		char b = '-';
		char c = '-';
		char d = '-';
		char e = '-';
		char f = '-';
		char g = '-';
		char h = '-';
		char i = '-';
		char j = '-';
	
		
	AminoAcidColumn col = new AminoAcidColumn(a, b, c, d, e, f, g, h, i, j);
	
	int nr = col.getOccurance('D');
	Assert.assertEquals(nr, 9);
	int res = col.howManyResidueTypes();
	Assert.assertEquals(res, 1);
	Assert.assertEquals(col.allButOneGaps(), true);
	
	
	
		
	}

}
