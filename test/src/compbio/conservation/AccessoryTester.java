package compbio.conservation;
import org.testng.Assert;
import org.testng.annotations.Test;

public class AccessoryTester {
	
	//@Test
	//public void  VoronoiTester() {
		
		//char a = 'D';
		//char b = 'D';
		//char c = 'D';
		//char d = 'D';
		//char e = 'D';
		//char f = 'E';
		///char g = 'E';
		//char h = 'E';
		//char i = 'F';
		//char j = 'F';
	
		//AminoAcidMatrix matrix = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i, j);
		
		//double result = ConservationAccessory.voronoiWeights(matrix, 1000)[0];
		
		//Assert.assertEquals(result, 1.0, 0.1);
	//}
	
	
	@Test
	public void normalizationTester() {
		
		double[] scores = {15, 13, 13, -32, 0};
		
		double[] normalized = ConservationAccessory.normalize01(scores);
		
		System.out.println("normal");
		
		for (int i = 0; i < normalized.length; i++) {
			
			System.out.println(normalized[i]);
		}
	}
	
	@Test
	public void inversedNormalizationTester() {
		
		double[] scores = {15, 13, 13, -32, 0};
		
		double[] normalized = ConservationAccessory.inversedNormalize01(scores);
		
		System.out.println("inversed");
		
		for (int i = 0; i < normalized.length; i++) {
			
			System.out.println(normalized[i]);
		}
	}

}
