package compbio.conservation;

import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.fail;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;
import compbio.data.sequence.UnknownFileFormatException;
import compbio.util.NullOutputStream;
import compbio.util.Timer;

public class ConservationTester {

	static ExecutorFactory efactory;

	@BeforeClass
	public void init() {
		ExecutorFactory.initExecutor(0,
				new PrintWriter(new NullOutputStream()),
				ExecutorFactory.ExecutorType.AsynchQueue);
	}

	@Test
	public void testMethods() {
		File input = new File(SlowMethodTester.DATA_PATH + File.separator
				+ SlowMethodTester.SMALL_AL);
		try {

			Conservation c = Conservation.getConservation(input, true,
					ExecutorFactory.getExecutor());
			Map<Method, double[]> results = c.calculateScores(EnumSet
					.allOf(Method.class));
			assertNotNull(results);
			c.printResults();
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

	@Test
	public void testSadler() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(
							SlowMethodTester.DATA_PATH + File.separator
									+ SlowMethodTester.SMALL_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(Method.SANDER);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			// ConservationScore2Tester.printScores(result, "Sander");
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}
}
