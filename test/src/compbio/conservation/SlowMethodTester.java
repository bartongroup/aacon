package compbio.conservation;

import static org.testng.Assert.fail;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.testng.annotations.Test;

import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;
import compbio.util.Timer;

public class SlowMethodTester {

	public final static String DATA_PATH = "test/data";

	static final String TINY_AL = "small.align";
	static final String SMALL_AL = "TO1296.fasta.align";
	static final String AVG_AL = "avg.aln.fa";
	static final String LARGE_AL = "1000x3000Dna.aln.fa";

	@Test(enabled = false)
	public void testSadler() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true);
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(Method.SANDER);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			ConservationScore2Tester.printScores(result, "Sander");
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test
	public void testSMERFS() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			// 675706 ms
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, false);
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.getSMERFS(7, SMERFSColumnScore.MID_SCORE,
					0.1);
			System.out.println("Calculating SMERFS scores: "
					+ timer.getStepTime());

			// this is a wrong call!
			System.out.println(Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	public static void main(String[] args) {
		/*
		 * int[][] a = new int[20000][5000]; System.out.println("T " +
		 * Runtime.getRuntime().totalMemory()); System.out.println("M " +
		 * Runtime.getRuntime().maxMemory()); System.out.println("F " +
		 * Runtime.getRuntime().freeMemory());
		 */
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			// 675706 ms
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, false);
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.getSMERFS(7, SMERFSColumnScore.MID_SCORE,
					0.1);
			System.out.println("Calculating SMERFS scores: "
					+ timer.getStepTime());

			// this is a wrong call!
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
