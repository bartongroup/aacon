package compbio.conservation;

import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

final class ExecutorFactory {

	public enum ExecutorType {
		SynchroniousCallerRuns, AsynchQueue
	};

	private static volatile ExecutorService executor;

	private static int threadNum;

	private static void initializeExecutor(int procNum, PrintWriter statWriter,
			ExecutorType etype) {
		int corenum = Runtime.getRuntime().availableProcessors();
		if (procNum < 1) {
			// Default - no thread number was set
			procNum = corenum;
		} else if (procNum > corenum * 2) {
			// To many cpus are defined -> user mistake
			String message = "Number of processors must be more than 1 and "
					+ "\n" + "less than the number of cores*2 " + "\n"
					+ "However given number of processors is " + procNum + "\n"
					+ "Changing number of processors to " + corenum
					+ " - the number of cores\n";
			System.err.println(message);
			statWriter.println(message);
			procNum = corenum;
		}
		statWriter.println("Using " + procNum + " CPUs");
		ExecutorFactory.threadNum = procNum;
		executor = getExecutor(etype);
	}

	static void initExecutor(int procNum, PrintWriter statWriter,
			ExecutorType etype) {
		if (executor == null) {
			synchronized (ExecutorFactory.class) {
				if (executor == null) {
					initializeExecutor(procNum, statWriter, etype);
				}
			}
		}
	}

	static ExecutorService getExecutor() {
		if (executor == null) {
			throw new IllegalStateException(
					"Please initialize the executor first!");
		}
		return executor;
	}

	private static ExecutorService getExecutor(ExecutorType etype) {
		if (executor == null) {
			switch (etype) {
			case AsynchQueue:
				executor = Executors.newFixedThreadPool(threadNum);
				break;
			case SynchroniousCallerRuns:
				executor = new ThreadPoolExecutor(threadNum, threadNum, 0L,
						TimeUnit.MILLISECONDS,
						new SynchronousQueue<Runnable>(),
						new ThreadPoolExecutor.CallerRunsPolicy());
				break;
			default:
				throw new RuntimeException("Unsupported executor type!");
			}
		}
		return executor;
	}

}
