package compbio.util;

import static org.testng.Assert.assertTrue;

import org.testng.annotations.Test;

public class UtilTester {

	@Test
	public void testRandomRange() {
		for (int i = 0; i < 10000; i++) {
			int num = Util.getRandomNumber(10, 99);
			assertTrue(num >= 10 && num <= 99);
		}
	}
}
