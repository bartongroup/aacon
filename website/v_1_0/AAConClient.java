/* Copyright (c) 2010 Peter Troshin
 *  
 *  Amino Acid Conservation Web Service client @version: 1.0     
 * 
 *  This library is free software; you can redistribute it and/or modify it under the terms of the
 *  Apache License version 2 as published by the Apache Software Foundation
 * 
 *  This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 *  even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache 
 *  License for more details.
 * 
 *  A copy of the license is in apache_license.txt. It is also available here:
 * @see: http://www.apache.org/licenses/LICENSE-2.0.txt
 * 
 * Any republication or derived work distributed in source code form
 * must include this copyright and license notice.
 */
package compbio.ws.client;

import static compbio.ws.client.Constraints.inputkey;
import static compbio.ws.client.Constraints.limitList;
import static compbio.ws.client.Constraints.outputkey;
import static compbio.ws.client.Constraints.paramFile;
import static compbio.ws.client.Constraints.paramList;
import static compbio.ws.client.Constraints.presetList;
import static compbio.ws.client.Constraints.presetkey;
import static compbio.ws.client.Constraints.pseparator;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.namespace.QName;
import javax.xml.ws.Service;
import javax.xml.ws.WebServiceException;

import compbio.data.msa.Annotation;
import compbio.data.sequence.FastaSequence;
import compbio.data.sequence.Score;
import compbio.data.sequence.SequenceUtil;
import compbio.data.sequence.UnknownFileFormatException;
import compbio.metadata.JobSubmissionException;
import compbio.metadata.Option;
import compbio.metadata.Preset;
import compbio.metadata.ResultNotAvailableException;
import compbio.metadata.WrongParameterException;
import compbio.ws.server.AAConWS;

/**
 * A command line client for AACon web service
 * 
 * @author pvtroshin
 * @version 1.0
 */
public class AAConClient {

	/**
	 * Web service host
	 */
	static final String hostname = "http://www.compbio.dundee.ac.uk/aacon";
	// static final String hostname = "http://localhost:8080/jabaws";

	/*
	 * Use java.util.Logger instead of log4j logger to reduce the size of the
	 * client package
	 */
	private static final Logger log = Logger.getLogger(AAConClient.class
			.getCanonicalName());

	/**
	 * The fully qualified web service namespace
	 */
	static final String QUALIFIED_SERVICE_NAME = "http://msa.data.compbio/01/12/2010/";

	/**
	 * Web service name
	 */
	static final String serviceName = "AAConWS";

	/**
	 * Calculate conservation for sequences loaded from the file
	 * 
	 * @param wsproxy
	 *            a web service proxy
	 * @param file
	 *            the file to read the results from
	 * @param preset
	 *            Preset to use optional
	 * @param customOptions
	 *            the list of options
	 * @return Set<Score> the conservation scores
	 * @throws UnknownFileFormatException
	 */
	static Set<Score> analize(File file, Annotation<AAConWS> wsproxy,
			Preset<AAConWS> preset, List<Option<AAConWS>> customOptions) {

		List<FastaSequence> fastalist = null;
		Set<Score> scores = null;
		try {
			fastalist = SequenceUtil.openInputStream(file.getAbsolutePath());

			String jobId = null;
			if (customOptions != null && preset != null) {
				System.out
						.println("WARN: Parameters (-f) are defined together with a preset (-r) ignoring preset!");
			}
			if (customOptions != null) {
				jobId = wsproxy.customAnalize(fastalist, customOptions);
			} else if (preset != null) {
				jobId = wsproxy.presetAnalize(fastalist, preset);
			} else {
				jobId = wsproxy.analize(fastalist);
			}
			Thread.sleep(1000);
			scores = wsproxy.getAnnotation(jobId);

		} catch (IOException e) {
			System.err
					.println("Exception while reading the input file. "
							+ "Check that the input file contains a list of fasta formatted sequences! "
							+ "Exception details are below:");
			e.printStackTrace();
		} catch (JobSubmissionException e) {
			System.err
					.println("Exception while submitting job to a web server. "
							+ "Exception details are below:");
			e.printStackTrace();
		} catch (ResultNotAvailableException e) {
			System.err.println("Exception while waiting for results. "
					+ "Exception details are below:");
			e.printStackTrace();
		} catch (InterruptedException ignored) {
			// ignore and propagate an interruption
			Thread.currentThread().interrupt();
		} catch (WrongParameterException e) {
			System.err
					.println("Exception while parsing the web method input parameters. "
							+ "Exception details are below:");
			e.printStackTrace();
		} catch (UnknownFileFormatException e) {
			System.err
					.println("Exception while attempting to read the input file "
							+ "Exception details are below:");
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
		return scores;
	}

	/**
	 * Connects to a AACon web service by the host and the service name
	 * 
	 * 
	 * @return {@link Annotation}
	 * @throws WebServiceException
	 *             if cannot connect to a web service
	 */
	public static Annotation<AAConWS> connect() throws WebServiceException {
		URL url = null;
		log.log(Level.FINE, "Attempting to connect...");
		try {
			url = new URL(hostname + "/" + "AAConWS" + "?wsdl");
		} catch (MalformedURLException e) {
			e.printStackTrace();
			// ignore as the host name is already verified
		}
		QName qname = new QName(QUALIFIED_SERVICE_NAME, "AAConWS");
		Service serv = Service.create(url, qname);
		QName portName = new QName(QUALIFIED_SERVICE_NAME, "AAConWS" + "Port");

		@SuppressWarnings("unchecked")
		Annotation<AAConWS> serviceIF = serv
				.getPort(portName, Annotation.class);

		log.log(Level.FINE, "Connected successfully!");
		return serviceIF;
	}

	/**
	 * Starts command line client, if no parameters are supplied prints help.
	 * 
	 * @param args
	 *            Usage: <Class or Jar file name> ACTION [OPTIONS]
	 * 
	 *            -i=<inputFile> - full path to fasta or Clustal formatted
	 *            alignment file
	 * 
	 *            -parameters - lists parameters supported by web service
	 * 
	 *            -presets - lists presets supported by web service
	 * 
	 *            -limits - lists web services limits. Please note that if input
	 *            file is specified other actions are ignored
	 * 
	 *            OPTIONS: (only for use with -i action):
	 * 
	 *            -r=<presetName> - name of the preset to use
	 * 
	 *            -o=<outputFile> - full path to the file where to write results
	 *            -f=<parameterInputFile> - the name of the file with the list
	 *            of parameters to use. Please note that -r and -f options
	 *            cannot be used together. Conservation is calculated with
	 *            either a preset or parameters from the file, but not both!
	 * 
	 */
	public static void main(String[] args) {

		if (args == null) {
			printUsage(1);
		}
		if (args.length < 1) {
			System.out.println("No options is specified! ");
			printUsage(1);
		}

		try {
			new AAConClient(args);
		} catch (IOException e) {
			log.log(Level.SEVERE, "IOException in client! " + e.getMessage(),
					e.getCause());
			System.err.println("Cannot write output file! Stack trace: ");
			e.printStackTrace();
		}
	}

	/**
	 * Prints AAConClient usage information to standard out
	 * 
	 * @param exitStatus
	 */
	static void printUsage(int exitStatus) {
		System.out.println();
		System.out.println("Usage: <Class or Jar file name> "
				+ " ACTION [OPTIONS] ");
		System.out.println();
		System.out.println("ACTIONS: ");
		System.out
				.println(inputkey
						+ pseparator
						+ "<inputFile> - full path to fasta or Clustal formatted alignment file ");
		System.out.println(paramList
				+ " - lists parameters supported by web service");
		System.out.println(presetList
				+ " - lists presets supported by web service");
		System.out.println(limitList + " - lists web services limits");
		System.out
				.println("Please note that if input file is specified other actions are ignored");

		System.out.println();
		System.out.println("OPTIONS (only for use with -i action):");

		System.out.println(presetkey + pseparator
				+ "<presetName> - name of the preset to use");
		System.out
				.println(outputkey
						+ pseparator
						+ "<outputFile> - full path to the file where to write the result");
		System.out
				.println("-f=<parameterInputFile> - the name of the file with the list of parameters to use.");
		System.out
				.println("Please note that -r and -f options cannot be used together. "
						+ "Conservation is calculated with either a preset or "
						+ "the parameters from the file, but not both!");

		System.exit(exitStatus);
	}

	/**
	 * Outputs AAcon results into the file represented by the outStream
	 * 
	 * @param outStream
	 * @param result
	 *            the AACon scores to output
	 */
	static void writeOut(OutputStream outStream, Set<Score> result) {
		try {
			Score.write(result, outStream);
		} catch (IOException e) {
			System.err
					.println("Problems writing output file! Stack trace is below: ");
			e.printStackTrace();
		} finally {
			if (outStream != null) {
				try {
					outStream.close();
				} catch (IOException ignored) {
					// e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Connects to the service and do the job as requested, if something goes
	 * wrong reports or/and prints usage help.
	 * 
	 * @param cmd
	 *            command line options
	 * @throws IOException
	 *             if the system cannot read/write from/into the file system
	 */
	@SuppressWarnings("unchecked")
	AAConClient(String[] cmd) throws IOException {

		File inputFile = IOHelper.getFile(cmd, inputkey, true);
		File outFile = IOHelper.getFile(cmd, outputkey, false);
		File parametersFile = IOHelper.getFile(cmd, paramFile, true);
		String presetName = CmdHelper.getPresetName(cmd);

		Annotation<AAConWS> msaws = connect();
		Preset<AAConWS> preset = null;
		if (presetName != null) {
			preset = MetadataHelper.getPreset(msaws, presetName);
		}
		List<Option<AAConWS>> customOptions = null;
		if (parametersFile != null) {
			List<String> prms = IOHelper.loadParameters(parametersFile);
			customOptions = MetadataHelper.processParameters(prms,
					msaws.getRunnerOptions());
		}
		Set<Score> result = null;
		if (inputFile != null) {
			System.out.println("Calculating conservation...");
			result = analize(inputFile, msaws, preset, customOptions);
			OutputStream outStream = null;
			if (outFile != null) {
				outStream = IOHelper.getOutStream(outFile);
			} else {
				// this stream is going to be closed later which is fine as
				// std.out will not be
				outStream = System.out;
			}
			writeOut(outStream, result);
			// stream is closed in the method no need to close it here
		}

		boolean listParameters = CmdHelper.listParameters(cmd);
		if (listParameters) {
			System.out.println(MetadataHelper.getParametersList(msaws));
		}
		boolean listPreset = CmdHelper.listPresets(cmd);
		if (listPreset) {
			System.out.println(MetadataHelper.getPresetList(msaws));
		}
		boolean listLimits = CmdHelper.listLimits(cmd);
		if (listLimits) {
			System.out.println(MetadataHelper.getLimits(msaws));
		}
		log.fine("Disconnecting...");
		((Closeable) msaws).close();
		log.fine("Disconnected successfully!");
	}
}
