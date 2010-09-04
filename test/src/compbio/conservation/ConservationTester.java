package compbio.conservation;

import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.fail;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;
import compbio.data.sequence.UnknownFileFormatException;

public class ConservationTester {

	static ExecutorFactory efactory;
	private Map<Method, double[]> norm_results = null;
	private Map<Method, double[]> results = null;

	private static File input = new File(SlowMethodTester.DATA_PATH
			+ File.separator + SlowMethodTester.AVG_AL4);

	@BeforeClass
	public void init() {
		// Make serial calculation
		ExecutorFactory.initExecutor(1);

		List<FastaSequence> sequences = null;
		try {
			sequences = SequenceUtil.readFasta(new FileInputStream(input));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
		assertNotNull(sequences);
		AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
		Conservation scores = new Conservation(alignment, true,
				ExecutorFactory.getExecutor());
		EnumSet<Method> set = EnumSet.allOf(Method.class); // EnumSet<Method>
		// set = EnumSet.range(Method.JORES, Method.SANDER);
		norm_results = scores.calculateScores(set);
		scores = new Conservation(alignment, false,
				ExecutorFactory.getExecutor());
		results = scores.calculateScores(set);
		// Shutdown the executor and complete submitted tasks
		ExecutorFactory.getExecutor().shutdown();
		try {
			ExecutorFactory.getExecutor().awaitTermination(5, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	@Test()
	public void testSingleMethods() {
		// Re-initiate executor with multiple threads
		ExecutorFactory.initExecutor(0);
		Conservation cons;
		try {
			cons = Conservation.getConservation(input, false,
					ExecutorFactory.getExecutor());
			double[] result = cons.calculateScore(Method.VALDAR);

			cons = Conservation.getConservation(input, false,
					ExecutorFactory.getExecutor());
			double[] result2 = cons.calculateScore(Method.VALDAR);

			Assert.assertEquals(results.get(Method.VALDAR).length,
					result.length);

			System.out.println(Arrays.toString(result));
			// System.out.println(Arrays.toString(result2));
			// System.out.println(Arrays.toString(results.get(Method.VALDAR)));

			Assert.assertTrue(Arrays.equals(result2, result));

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (UnknownFileFormatException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

	@Test()
	public void testMethods() {
		try {

			Conservation c = Conservation.getConservation(input, false,
					ExecutorFactory.getExecutor());
			Map<Method, double[]> apiresults = c.calculateScores(EnumSet
					.allOf(Method.class));
			assertNotNull(results);
			for (Method method : apiresults.keySet()) {
				if (method == Method.LANDGRAF || method == Method.VALDAR
						|| method == Method.KABAT) {
					// Landgrap results never repeats as they have random
					// element
					continue;
				}
				double[] result = results.get(method);
				double[] apiresult = apiresults.get(method);
				Assert.assertNotNull(apiresult);
				System.out.println(Arrays.toString(result));
				System.out.println(Arrays.toString(apiresult));

				Assert.assertTrue(Arrays.equals(apiresult, result),
						"Methods results: " + method.toString()
								+ " is not equal!");
			}
			c.printResults(Format.RESULT_NO_ALIGNMENT);
			// System.out.println(results);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getMessage());
		} catch (UnknownFileFormatException e) {
			e.printStackTrace();
			fail(e.getMessage());
		}
	}

}
