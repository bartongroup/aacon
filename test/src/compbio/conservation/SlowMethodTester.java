package compbio.conservation;

import static org.testng.Assert.fail;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.SequenceUtil;
import compbio.util.NullOutputStream;
import compbio.util.Timer;

public class SlowMethodTester {

	public final static String DATA_PATH = "test/data";

	static final String TINY_AL = "small.align";
	static final String SMALL_AL = "TO1296.fasta.align";
	static final String AVG_AL = "avg.aln.fa";
	static final String LARGE_AL = "1000x3000Dna.aln.fa";
	static ExecutorFactory efactory;

	@BeforeClass
	public void init() {
		ExecutorFactory.initExecutor(0,
				new PrintWriter(new NullOutputStream()),
				ExecutorFactory.ExecutorType.AsynchQueue);
	}

	@Test()
	public void testSadler() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, false,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(Method.SANDER);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			// Conservation.printResults(result, Method.SANDER);
			scores.outputResults(new File("results.txt"),
					Format.RESULT_NO_ALIGNMENT);
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test(enabled = false)
	public void testLandgraf() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(Method.LANDGRAF);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#LADGRAF " + Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test(enabled = false)
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

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.getSMERFS(11, SMERFSColumnScore.MID_SCORE,
					0.2);
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

	@Test()
	public void reproduceIssue1() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(Method.JORES);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#KABAT " + Arrays.toString(result));
			System.out.println("Total: " + timer.getTotalTime());

		} catch (FileNotFoundException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		} catch (IOException e) {
			e.printStackTrace();
			fail(e.getLocalizedMessage());
		}
	}

	@Test()
	public void reproduceIssue4() {
		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(Method.JORES);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#JORES " + Arrays.toString(result));
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
		ExecutorFactory.initExecutor(0,
				new PrintWriter(new NullOutputStream()),
				ExecutorFactory.ExecutorType.AsynchQueue);

		try {
			Timer timer = new Timer(TimeUnit.MILLISECONDS);
			List<FastaSequence> sequences = SequenceUtil
					.readFasta(new FileInputStream(new File(DATA_PATH
							+ File.separator + AVG_AL)));
			System.out.println("Loading sequences: " + timer.getStepTime());

			AminoAcidMatrix alignment = new AminoAcidMatrix(sequences, null);
			System.out.println("Converting to Matrix: " + timer.getStepTime());

			Conservation scores = new Conservation(alignment, true,
					ExecutorFactory.getExecutor());
			System.out.println("Constructing conservation scores: "
					+ timer.getStepTime());

			double[] result = scores.calculateScore(Method.JORES);
			System.out.println("Calculating sadler scores: "
					+ timer.getStepTime());

			System.out.println("#KABAT " + Arrays.toString(result));
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
