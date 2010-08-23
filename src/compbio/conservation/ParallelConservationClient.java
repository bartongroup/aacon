package compbio.conservation;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import compbio.data.sequence.FastaSequence;
import compbio.util.Timer;

/**
 * TODO to complete
 * 
 * add n-cpu parameter
 * 
 * @author pvtroshin
 */
public final class ParallelConservationClient {

	private final Map<Method, double[]> results = new EnumMap<Method, double[]>(
			Method.class);
	private static volatile ExecutorService executor;

	private void initExecutor(int procNum) {

		int corenum = Runtime.getRuntime().availableProcessors();
		if (procNum <= 1 || procNum > corenum * 2) {
			System.err.println("Number of processors must be more than 1 and "
					+ "less than the number of cores*2"
					+ "However given number of processors is " + procNum);
			System.err.println("Changing number of processors to " + corenum
					+ " - the number of cores");
			procNum = corenum;
		}
		if (executor == null) {
			synchronized (ParallelConservationClient.class) {
				if (executor == null) {
					executor = new ThreadPoolExecutor(procNum, procNum, 0L,
							TimeUnit.MILLISECONDS,
							new SynchronousQueue<Runnable>(),
							new ThreadPoolExecutor.CallerRunsPolicy());
				}
			}
		}
	}

	public ParallelConservationClient(String[] cmd) throws IOException,
			InterruptedException {

		Timer timer = Timer.getMilliSecondsTimer();
		initExecutor(4);

		// default values
		int SMERFSWidth = SMERFSColumnScore.DEFAULT_WINDOW_SIZE;
		SMERFSColumnScore score = SMERFSColumnScore.MID_SCORE;
		double SMERFSGapTreshold = SMERFSColumnScore.DEFAULT_GAP_THRESHOLD;

		Set<Method> methods = CmdParser.getMethodNames(cmd);

		String inFilePath = CmdParser.getInputFilePath(cmd);

		if (methods != null && inFilePath != null) {
			String format = CmdParser.getFormat(cmd);
			String outFilePath = CmdParser.getOutputFilePath(cmd);

			Format outFormat = Format.getFormat(format);
			String[] SMERFSDetails = CmdParser.getSMERFSDetails(cmd);
			if (SMERFSDetails != null) {
				if (SMERFSDetails.length == 3) {
					try {
						SMERFSWidth = Integer.parseInt(SMERFSDetails[0]);
					} catch (NumberFormatException e) {
						SMERFSWidth = -1;
					}
					score = SMERFSColumnScore
							.getSMERFSColumnScore(SMERFSDetails[1]);
					try {
						SMERFSGapTreshold = Double
								.parseDouble(SMERFSDetails[2]);
					} catch (NumberFormatException e) {
						SMERFSGapTreshold = -1;
					}
				} else {
					System.err.println("To run SMERFS three arguments are"
							+ " needed, window width, how to give "
							+ "scores to columns and a gap treshold.");
					System.exit(1);
				}
			}
			boolean normalize = CmdParser.getNormalize(cmd);
			String[] gap = CmdParser.getGapChars(cmd);
			Character[] gapChars = CmdParser.extractGapChars(gap);

			String statFile = CmdParser.getStatFilePath(cmd);
			if (statFile == null) {
				timer.setStatOutput(null);
			} else {
				timer.setStatOutput(new FileOutputStream(statFile));
			}

			List<FastaSequence> sequences = CmdParser
					.openInputStream(inFilePath);

			if (sequences != null) {
				AminoAcidMatrix alignment = new AminoAcidMatrix(sequences,
						gapChars);
				timer.println("Start time: " + CmdParser.getDateTime());
				timer.println("Alignment loaded in: " + timer.getStepTime()
						+ " ms");
				timer.println("Alignment has: " + alignment.numberOfRows()
						+ " sequences.");

				ConservationScores2 scores = new ConservationScores2(alignment);

				MethodWrapper wrapper = null;
				List<MethodWrapper> tasks = new ArrayList<MethodWrapper>();

				for (Method method : methods) {
					wrapper = new MethodWrapper(method, scores, normalize,
							timer);
					if (method == Method.SMERFS && SMERFSDetails != null) {
						// this meant to override
						wrapper = new MethodWrapper(scores, normalize,
								SMERFSWidth, score, SMERFSGapTreshold, timer);
					}
					tasks.add(wrapper);
				}
				List<Future<MethodWrapper>> rawResults = executor
						.invokeAll(tasks);
				for (Future<MethodWrapper> rawResult : rawResults) {
					MethodWrapper entry = null;
					try {
						entry = rawResult.get();
					} catch (ExecutionException e) {
						// TODO better exception handling
						e.printStackTrace();
					}
					results.put(entry.method, entry.conservation);
				}
				executor.shutdown();

				ConservationFormatter.formatResults(results, outFilePath,
						outFormat, alignment);

				timer.println("Total calculation time: "
						+ timer.getTotalTime(TimeUnit.SECONDS) + " s");
				timer.println("End time: " + CmdParser.getDateTime());
				timer.getStatWriter().close();
			} else {
				System.out.println("No input found in " + inFilePath
						+ " ! Exiting");
			}
		}

	}

	public static void main(String[] args) {

		if (args == null) {
			System.out.println("No parameters were suppled");
			System.out.println();
			System.out.print(CmdParser.CONSERVATION_HELP);
		}
		if (args.length < 2) {
			System.out
					.println("Method names, input file paths are required. Application will"
							+ " not run until these 2 arguments are provided.");
			System.out
					.println("If you want results printed, both format an input "
							+ "file path have to be provided");
			System.out.println();
			System.out.print(CmdParser.CONSERVATION_HELP);
		}
		try {
			ParallelConservationClient cons = new ParallelConservationClient(
					args);
		} catch (IOException e) {
			System.err.println("Fail to write to the file system! "
					+ e.getLocalizedMessage());
			e.printStackTrace();
		} catch (InterruptedException e) {
			System.err.println("Interrupted!");
			e.printStackTrace();
		}
	}
}
