package compbio.conservation;

import org.testng.annotations.Test;
import java.util.*;

public class ConservationScoreTester {

	static void printScores(double[] results, String tag) {
		
		System.out.println(tag);
		
		for(int i = 0; i < results.length; i++ ) {
			
			System.out.println(results[i]);
		}
	}
	
	char a = 'D';
	char b = 'E';
	char c = 'P';
	
	char d = 'D';
	char e = '-';
	char f = 'E';
	
	char g = 'D';
	char h = '-';
	char i = 'K';
	
	AminoAcidMatrix alignment = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i);
	
	ConservationScores scores = new ConservationScores(alignment);
	
	public void kabatTester() {
		
		double[] result = scores.kabatScore(true);
		
		printScores(result, "Kabat");
		
	}
	
	@Test
	public void totalAcidsWillSetsTester() {
		
		Map<String, Integer> sets = alignment.totalAcidsWillSets();
		
		Iterator<String> itr = sets.keySet().iterator();
		
		while(itr.hasNext()) {
			
			String key = itr.next();
			
			System.out.println("Key" + key + "Value" + sets.get(key));
		}
		
		
	}
	
	@Test
	public void joresTester() {
		
		double[] result = scores.joresScore(true);
		
		printScores(result, "Jores");
		
	}
	
	@Test
	public void schneiderTester() {
		
		double[] result = scores.schneiderScore(true);
		
		printScores(result, "Schneider");
		
	}
	
	@Test
	public void shenkinTester() {
		
		double[] result = scores.shenkinScore(true);
		
		printScores(result, "Shenkin");
		
	}
	
	@Test
	public void gersteinTester() {
		
		double[] result = scores.gersteinScore(true);
		
		printScores(result, "Gerstein");
		
	}
	
	@Test
	public void smallestTaylosetGapsTester() {
		
		double[] result = scores.SmallestTaylorSetGaps(true);
		
		printScores(result, "SmallestTaylorSetGaps");
		
	}
	
	@Test
	public void smallestTaylorSetNoGapsTester() {
		
		double[] result = scores.SmallestTaylorSetNoGaps(true);
		
		printScores(result, "SmallestTaylorSetNoGaps");
		
	}
	
	@Test
	public void zvelibilTester() {
		
		double[] result = scores.zvelibilScore(true);
		
		printScores(result, "Zvelibil");
		
	}
	
	@Test
	public void karlinTester() {
		
		double[] result = scores.karlinScore(true);
		
		printScores(result, "Karlin");
		
	}
	
	@Test
	public void armonTester() {
		
		double[] result = scores.armonScore(true);
		
		printScores(result, "Armon");
		
	}
	
	@Test
	public void thompsonTester() {
		
		double[] result = scores.thompsonScore(true);
		
		printScores(result, "Thompson");
		
	}
	
	@Test
	public void notLancetTester() {
		
		double[] result = scores.notLancetScore(true);
		
		printScores(result, "NotLancet");
		
	}
	
	@Test
	public void mirnyTester() {
		
		double[] result = scores.mirnyScore(true);
		
		printScores(result, "Mirny");
		
	}
	
	@Test
	public void williamsonTester() {
		
		double[] result = scores.williamsonScore(true);
		
		printScores(result, "Williamson");
		
	}
	
	@Test
	public void landgrafTester() {
		
		double[] result = scores.landgrafScore(true);
		
		printScores(result, "Landgraf");
		
	}
	
	@Test
	public void sanderTester() {
		
		double[] result = scores.sanderScore(true);
		
		printScores(result, "Sander");
		
	}
	
	@Test
	public void valdarTester() {
		
		double[] result = scores.valdarScore(true);
		
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
