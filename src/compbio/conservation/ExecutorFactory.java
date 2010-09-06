/*
 * Copyright (c) 2010 Agnieszka Golicz & Peter Troshin 
 * 
 * Amino Acid Conservation @version: 1.0 
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the Apache License version 2 as published by the
 * Apache Software Foundation This library is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache
 * License for more details. A copy of the license is in apache_license.txt. It
 * is also available here: http://www.apache.org/licenses/LICENSE-2.0.txt 
 * Any republication or derived work distributed in source code form must 
 * include this copyright and license notice.
 * 
 */
package compbio.conservation;

import java.io.PrintWriter;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.SynchronousQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import compbio.util.NullOutputStream;

/**
 * 
 * Constructs and stores an instance of the {@link ExecutorService} for use by
 * the clients. For AA Conservation calculation, the best executor to use is the
 * one with the fixed thread number equals to the number of CPU cores and
 * Blocking queue as the task queue.
 * 
 * Executors.newFixedThreadPool(threadNum);
 * 
 * @author Peter Troshin
 * 
 */
public final class ExecutorFactory {

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
	 * Initializes the executor and stores its instance in the local static
	 * field. Only one instance of the executor is kept at any time.
	 * 
	 * If the executor is re-initialized with a different number of threads
	 * (procNum) then it is up to the caller to make sure that that all previous
	 * scheduled tasks have completed. Something like the following can be used
	 * to ensure that all tasks have been completed:
	 * 
	 * <pre>
	 * executor.shutdown();
	 * 
	 * executor.awaitTermination(60, TimeUnit.MINUTES);
	 * 
	 * </pre>
	 * 
	 * This code will cause the executor to wait for a maximum of 1 hour for all
	 * scheduled tasks to complete.
	 * 
	 * If this method call with the executor instance already initialized, then
	 * executor.shutdownNow() method will be invoked and all unfinished tasks
	 * are cancelled.
	 * 
	 * @param procNum
	 *            the number of thread to initialize executor with. For the
	 *            purpose of AA conservation calculations for the maximum
	 *            performance it is recommended that the {@code procNum} is
	 *            equal to the number of CPU cores. procNum must not exceed
	 *            2xcores. If it does, then it will be initialized with the n
	 *            cores equal to the number of CPU cores on the machine and the
	 *            error message printed to the standard error stream.
	 * @param statWriter
	 *            the PrintWriter to dump the execution statistics to.
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

	/**
	 * Initializes the instance of the {@link ExecutorService} with a number of
	 * computing thread equals to {@code procNum}. The execution statistics is
	 * voided.
	 * 
	 * @param procNum
	 *            the number of threads to use
	 */
	public static void initExecutor(int procNum) {
		initExecutor(procNum, new PrintWriter(new NullOutputStream()));
	}

	/**
	 * Initializes the instance of the {@link ExecutorService} with the number
	 * of threads equals to the number of cores available on the computer. The
	 * execution statistics is voided.
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
	 *            the PrintWriter to dump the execution statistics into
	 */
	public static void initExecutor(PrintWriter statWriter) {
		initExecutor(0, statWriter);
	}

	/**
	 * Return the executor which were initialized. Please note that one of the
	 * init methods of this class must be called to initialize the executor
	 * before this method!
	 * 
	 * @return the ExecutorService
	 * @throws IllegalStateException
	 *             if the executor has not been initialized yet
	 */
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
