package compbio.conservation;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.regex.Pattern;

import compbio.data.sequence.FastaSequence;

public class ResultReader {

	/**
	 * Reads printed results from the file.
	 * 
	 * @param inStream
	 *            input stream
	 * @return map mapping a map of results(with method names as keys and result
	 *         arrays as values) to the list of FASTA sequences that constitute
	 *         the alignment.
	 * @throws IOException
	 */
	static Map<Map<Method, double[]>, List<FastaSequence>> readFile(
			InputStream inStream) throws IOException {

		Map<Map<Method, double[]>, List<FastaSequence>> result = new HashMap<Map<Method, double[]>, List<FastaSequence>>();
		Map<Method, double[]> resultMap = new EnumMap<Method, double[]>(
				Method.class);
		BufferedReader inResults = new BufferedReader(new InputStreamReader(
				inStream));
		List<FastaSequence> seqList = new ArrayList<FastaSequence>();
		Pattern pattern = Pattern.compile("\\s+");
		String line;
		String lineString = null;
		// line = inResults.readLine();
		// line = inResults.readLine();
		do {
			line = inResults.readLine();
			// System.out.println(line);
			if (line == null || line.startsWith("#") || line.startsWith(">")) {
				if (lineString != null) {
					if (lineString.startsWith("#")) {
						parseResults(lineString, resultMap, pattern);
					}
					if (lineString.startsWith(">")) {
						parseSequences(lineString, seqList, pattern);
					}
				}
				lineString = line;
			} else {
				lineString += line;
			}
		} while (line != null);
		result.put(resultMap, seqList);
		return result;
	}

	static void parseResults(String resultStr, Map<Method, double[]> resultMap,
			Pattern pattern) {

		String resultStrTemp = pattern.matcher(resultStr.trim())
				.replaceAll(" ");
		String[] results = resultStrTemp.split(" ");
		String name = results[0].substring(1);
		System.out.println(name);
		double[] resultsNum = new double[results.length - 1];
		for (int i = 0; i < resultsNum.length; i++) {
			resultsNum[i] = Double.parseDouble(results[i + 1]);
			System.out.println(resultsNum[i]);
		}
		resultMap.put(Method.getMethod(name), resultsNum);
	}

	static void parseSequences(String lineStr, List<FastaSequence> list,
			Pattern pattern) {

		StringTokenizer tokens = new StringTokenizer(lineStr, " ");
		String name = tokens.nextToken().trim();
		System.out.println(name);
		Pattern pattern2 = Pattern.compile(name);
		String seqStrMod = pattern.matcher(
				pattern2.matcher(lineStr).replaceFirst("")).replaceAll("");
		System.out.println(seqStrMod);
		list.add(new FastaSequence(name, seqStrMod));
	}

}
