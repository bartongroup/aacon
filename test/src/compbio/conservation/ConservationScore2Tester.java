package compbio.conservation;

import org.testng.annotations.Test;
import java.util.*;

public class ConservationScore2Tester {

	static void printScores(double[] results, String tag) {
		
		System.out.println(tag);
		
		for(int i = 0; i < results.length; i++ ) {
			
			System.out.println(results[i]);
		}
	}
	
	char a = 'D';
	char b = 'F';
	char c = 'P';
	
	char d = 'D';
	char e = 'K';
	char f = 'P';
	
	char g = 'D';
	char h = '-';
	char i = 'K';
	
	
	
	AminoAcidMatrix alignment = new AminoAcidMatrix(a, b, c, d, e, f, g ,h,i);
	
	ConservationScores2 scores = new ConservationScores2(alignment);
	
	@Test
	public void kabatTester() {
		
		double[] result = scores.calculateScore(Method.KABAT, true);
		
		printScores(result, "Kabat");
		
	}
	
	@Test
	public void totalAcidsWillSetsTester() {
		
		Map<String, Integer> sets = alignment.totalAcidsWillSets();
		
		Iterator<String> itr = sets.keySet().iterator();
		
		while(itr.hasNext()) {
			
			String key = itr.next();
			
			//System.out.println("Key" + key + "Value" + sets.get(key));
		}
		
		
	}
	
	@Test
	public void joresTester() {
		
		double[] result = scores.calculateScore(Method.JORES, true);
		
		printScores(result, "Jores");
		
	}
	
	@Test
	public void schneiderTester() {
		
		double[] result = scores.calculateScore(Method.SCHNEIDER, true);
		
		printScores(result, "Schneider");
		
	}
	
	@Test
	public void shenkinTester() {
		
		double[] result = scores.calculateScore(Method.SHENKIN, true);
		
		printScores(result, "Shenkin");
		
	}
	
	@Test
	public void gersteinTester() {
		
		double[] result = scores.calculateScore(Method.GERSTEIN, true);
		
		printScores(result, "Gerstein");
		
	}
	
	@Test
	public void smallestTaylosetGapsTester() {
		
		double[] result = scores.calculateScore(Method.TAYLOR_GAPS, true);
		
		printScores(result, "SmallestTaylorSetGaps");
		
	}
	
	@Test
	public void smallestTaylorSetNoGapsTester() {
		
		double[] result = scores.calculateScore(Method.TAYLOR_NO_GAPS, true);
		
		printScores(result, "SmallestTaylorSetNoGaps");
		
	}
	
	@Test
	public void zvelibilTester() {
		
		double[] result = scores.calculateScore(Method.ZVELIBIL, true);
		
		printScores(result, "Zvelibil");
		
	}
	
	@Test
	public void karlinTester() {
		
		double[] result = scores.calculateScore(Method.KARLIN, true);
		
		printScores(result, "Karlin");
		
	}
	
	@Test
	public void armonTester() {
		
		double[] result = scores.calculateScore(Method.ARMON, true);
		
		printScores(result, "Armon");
		
	}
	
	@Test
	public void thompsonTester() {
		
		double[] result = scores.calculateScore(Method.THOMPSON, true);
		
		printScores(result, "Thompson");
		
	}
	
	@Test
	public void notLancetTester() {
		
		double[] result = scores.calculateScore(Method.NOT_LANCET, true);
		
		printScores(result, "NotLancet");
		
	}
	
	@Test
	public void mirnyTester() {
		
		double[] result = scores.calculateScore(Method.MIRNY, true);
		
		printScores(result, "Mirny");
		
	}
	
	@Test
	public void williamsonTester() {
		
		double[] result = scores.calculateScore(Method.WILLIAMSON, true);
		
		printScores(result, "Williamson");
		
	}
	
	@Test
	public void landgrafTester() {
		
		double[] result = scores.calculateScore(Method.LANDGRAF, true);
		
		printScores(result, "Landgraf");
		
	}
	
	@Test
	public void sanderTester() {
		
		double[] result = scores.calculateScore(Method.SANDER, true);
		
		printScores(result, "Sander");
		
	}
	
	@Test
	public void valdarTester() {
		
		double[] result = scores.calculateScore(Method.VALDAR, true);
		
		printScores(result, "Valdar");
		
	}
	
	@Test
	public void tester() {
		
		List<Map<Character, Integer>> map = alignment.getTotalAcidsFreqByCol();
		
		System.out.println("freqs");
		
		for (int i = 0; i < map.size(); i++) {
			
			Iterator<Character> itr = map.get(i).keySet().iterator();
			
			while(itr.hasNext()) {
				
				Character key = itr.next();
				
				System.out.println(i + " key " + key + " value " + map.get(i).get(key));
			}
			
			
		}
	}
}
