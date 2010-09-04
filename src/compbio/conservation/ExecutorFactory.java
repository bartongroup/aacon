package compbio.conservation;

import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import compbio.util.NullOutputStream;

final class ExecutorFactory {

	private enum ExecutorType {
		SynchroniousCallerRuns, AsynchQueue
	};

	private static volatile ExecutorService executor;

	private static volatile int threadNum = -1;

	private static void initializeExecutor(int procNum, PrintWriter statWriter) {
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
		executor = getExecutor(ExecutorType.AsynchQueue);
	}

	/**
	 * If the executor is re-initialized with a different number of threads
	 * (procNum) then it is up to the caller to make sure that that all previous
	 * scheduled tasks have completed. Upon this method call
	 * executor.shutdownNow() method will be invoked and all unfinished tasks
	 * will be cancelled.
	 * 
	 * @param procNum
	 * @param statWriter
	 */
	public static void initExecutor(int procNum, PrintWriter statWriter) {
		if (executor == null || threadNum != procNum) {
			synchronized (ExecutorFactory.class) {
				if (threadNum != procNum) {
					if (executor != null) {
						// make sure all scheduled tasks are completed
						executor.shutdownNow();
					}
					initializeExecutor(procNum, statWriter);
				}
			}
		}
	}

	public static void initExecutor(int procNum) {
		initExecutor(procNum, new PrintWriter(new NullOutputStream()));
	}

	/**
	 * Initializes the executor with the number of threads equals to the number
	 * of cores available on the computer. The execution statistics is voided.
	 */
	public static void initExecutor() {
		initExecutor(0, new PrintWriter(new NullOutputStream()));
	}

	/**
	 * Initializes the executor with the number of threads equals to the number
	 * of cores available on the computer. The execution statistics is dumped
	 * into the statWriter
	 * 
	 * @param statWriter
	 */
	public static void initExecutor(PrintWriter statWriter) {
		initExecutor(0, statWriter);
	}

	public static ExecutorService getExecutor() {
		if (executor == null) {
			throw new IllegalStateException(
					"Please initialize the executor first!");
		}
		return executor;
	}

	private static ExecutorService getExecutor(ExecutorType etype) {
		switch (etype) {
		case AsynchQueue:
			executor = Executors.newFixedThreadPool(threadNum);
			break;
		case SynchroniousCallerRuns:
			executor = new ThreadPoolExecutor(threadNum, threadNum, 0L,
					TimeUnit.MILLISECONDS, new SynchronousQueue<Runnable>(),
					new ThreadPoolExecutor.CallerRunsPolicy());
			break;
		default:
			throw new RuntimeException("Unsupported executor type!");
		}
		return executor;
	}

}
