package compbio.conservation;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.testng.annotations.Test;

public class ConservationScore2Tester {

	static void printScores(double[] results, String tag) {

		System.out.println(tag);

		for (int i = 0; i < results.length; i++) {

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

	AminoAcidMatrix alignment = new AminoAcidMatrix(a, b, c, d, e, f, g, h, i);

	Conservation scores = new Conservation(alignment, true);

	@Test
	public void kabatTester() {

		double[] result = scores.calculateScore(Method.KABAT);

		printScores(result, "Kabat");

	}

	@Test
	public void totalAcidsWillSetsTester() {

		Map<String, Integer> sets = alignment.totalAcidsWillSets();

		Iterator<String> itr = sets.keySet().iterator();

		while (itr.hasNext()) {

			String key = itr.next();

			// System.out.println("Key" + key + "Value" + sets.get(key));
		}

	}

	@Test
	public void joresTester() {

		double[] result = scores.calculateScore(Method.JORES);

		printScores(result, "Jores");

	}

	@Test
	public void schneiderTester() {

		double[] result = scores.calculateScore(Method.SCHNEIDER);

		printScores(result, "Schneider");

	}

	@Test
	public void shenkinTester() {

		double[] result = scores.calculateScore(Method.SHENKIN);

		printScores(result, "Shenkin");

	}

	@Test
	public void gersteinTester() {

		double[] result = scores.calculateScore(Method.GERSTEIN);

		printScores(result, "Gerstein");

	}

	@Test
	public void smallestTaylosetGapsTester() {

		double[] result = scores.calculateScore(Method.TAYLOR_GAPS);

		printScores(result, "SmallestTaylorSetGaps");

	}

	@Test
	public void smallestTaylorSetNoGapsTester() {

		double[] result = scores.calculateScore(Method.TAYLOR_NO_GAPS);

		printScores(result, "SmallestTaylorSetNoGaps");

	}

	@Test
	public void zvelibilTester() {

		double[] result = scores.calculateScore(Method.ZVELIBIL);

		printScores(result, "Zvelibil");

	}

	@Test
	public void karlinTester() {

		double[] result = scores.calculateScore(Method.KARLIN);

		printScores(result, "Karlin");

	}

	@Test
	public void armonTester() {

		double[] result = scores.calculateScore(Method.ARMON);

		printScores(result, "Armon");

	}

	@Test
	public void thompsonTester() {

		double[] result = scores.calculateScore(Method.THOMPSON);

		printScores(result, "Thompson");

	}

	@Test
	public void notLancetTester() {

		double[] result = scores.calculateScore(Method.NOT_LANCET);

		printScores(result, "NotLancet");

	}

	@Test
	public void mirnyTester() {

		double[] result = scores.calculateScore(Method.MIRNY);

		printScores(result, "Mirny");

	}

	@Test
	public void williamsonTester() {

		double[] result = scores.calculateScore(Method.WILLIAMSON);

		printScores(result, "Williamson");

	}

	@Test
	public void landgrafTester() {

		double[] result = scores.calculateScore(Method.LANDGRAF);

		printScores(result, "Landgraf");

	}

	@Test
	public void sanderTester() {

		double[] result = scores.calculateScore(Method.SANDER);

		printScores(result, "Sander");

	}

	@Test
	public void valdarTester() {

		double[] result = scores.calculateScore(Method.VALDAR);

		printScores(result, "Valdar");

	}

	@Test
	public void tester() {

		List<Map<Character, Integer>> map = alignment.getTotalAcidsFreqByCol();

		System.out.println("freqs");

		for (int i = 0; i < map.size(); i++) {

			Iterator<Character> itr = map.get(i).keySet().iterator();

			while (itr.hasNext()) {

				Character key = itr.next();

				System.out.println(i + " key " + key + " value "
						+ map.get(i).get(key));
			}

		}
	}
}
