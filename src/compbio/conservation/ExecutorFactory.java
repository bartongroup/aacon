package compbio.conservation;

import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

final class ExecutorFactory {

	private static volatile ExecutorService syncExecutor;
	private static volatile ExecutorService bQueueExecutor;

	private static volatile ExecutorFactory INSTANCE;

	private static int threadNum;

	private ExecutorFactory(int procNum, PrintWriter statWriter) {
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
	}

	static ExecutorFactory getFactory(int procNum, PrintWriter statWriter) {
		if (INSTANCE == null) {
			synchronized (ExecutorFactory.class) {
				if (INSTANCE == null) {
					INSTANCE = new ExecutorFactory(procNum, statWriter);
				}
			}
		}
		return INSTANCE;
	}

	ExecutorService getSynchroneousCallerRunsExecutor() {
		if (syncExecutor == null) {
			initSynchroneousCallerRunsExecutor();
		}
		return syncExecutor;
	}

	ExecutorService getQueueExecutor() {
		if (bQueueExecutor == null) {
			initBlockingExecutor();
		}
		return bQueueExecutor;
	}

	private static void initSynchroneousCallerRunsExecutor() {
		if (syncExecutor == null) {
			synchronized (ExecutorFactory.class) {
				if (syncExecutor == null) {
					if (threadNum == -1) {
						throw new IllegalStateException(
								"Must initialize the factory by calling initExecutorFactory() first!");
					}
					syncExecutor = new ThreadPoolExecutor(threadNum, threadNum,
							0L, TimeUnit.MILLISECONDS,
							new SynchronousQueue<Runnable>(),
							new ThreadPoolExecutor.CallerRunsPolicy());
				}
			}
		}
	}

	private static void initBlockingExecutor() {
		if (bQueueExecutor == null) {
			synchronized (ExecutorFactory.class) {
				if (bQueueExecutor == null) {
					bQueueExecutor = Executors.newFixedThreadPool(threadNum);
				}
			}
		}
	}

	void shutdownExecutors() {
		syncExecutor.shutdown();
		bQueueExecutor.shutdown();
	}
}
