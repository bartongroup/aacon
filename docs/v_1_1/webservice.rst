Web Service
===========

This page provides some tips on using AACon web service. Please note that for now all the examples are in `Java`_. Other languages will follow given a sufficient demand. Please also note that AACon web service client requires `Java`_ 6 or later.

------------

.. _web_details:

Web Service Details
-------------------

+------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Feature                | Description                                                                                                                                                                            |
+========================+========================================================================================================================================================================================+
| Operated by            | The University of Dundee, Computational Biology Group headed by Prof Geoff Barton. It is backed up by the College of Life Sciences HPC cluster.                                        |
+------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Service web site       | http://www.compbio.dundee.ac.uk/aacon                                                                                                                                                  |
+------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Execution limits       | Tasks exceeding 5000×1000 (sequences per letters) will not be accepted for alignment. If you would like to work with bigger alignments consider using a command line version of AACon. |
+------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| AACon web service WSDL | http://www.compbio.dundee.ac.uk/aacon/AAConWS?wsdl                                                                                                                                     |
+------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


.. tip:: For a more detailed description of all available types and their functions please refer to the `data model javadoc`_.


------------

.. _web_usage:

Using Web Service CLI Client
----------------------------

The command client can be used to calculate conservation using public AACon web service. The client is operating system independent and supports most of the AACon web service functions. Using this client you could calculate conservation using web service features such as presets and parameters. Below is the list of options supported by the command line client.

::

    Usage: java -jar <path_to_jar_file> ACTION [OPTIONS]
    ACTIONS:
    -i=<inputFile> - full path to Fasta or Clustal formatted alignment file
    -parameters - lists parameters supported by web service
    -presets - lists presets supported by web service
    -limits - lists web services limits
    Please note that if input file is specified other actions are ignored
    OPTIONS: (only for use with -i action):
    -r=<presetName> - name of the preset to use
    -o=<outputFile> - full path to the file where to write the results
    -f=<parameterInputFile> - the name of the file with the list of parameters to use.
    Please note that -r and -f options cannot be used together. Conservation can be calculated with either a preset or the parameters from the file, but not both!


Using the command line web service client is easy. For example to calculate conservation for the alignment from input.fasta file using AACon web service with default settings and print the results to the console. By default, the conservation is calculated by *SHENKIN* method.

.. code-block:: bash

    java -jar aacon-client.jar -i=d:\input.fasta


Content of input.fasta file is show below (please note sequences has been trimmed for clarity)

::

    >QUERY
    MEWLKIALVAFSGRGNDVDYPCSILGRVANNDKELVNIL
    >Q4VBK1
    ---------------PDTDEPVWILGACYNVKTKKSELL
    >A7RTH8
    ---MDAACVTYEGVSHETEEDVWILGKRYNIDMGYLNTD


Calculate conservation using 12 fast conservation methods (using preset "Quick conservation"), write output results to the result.txt file

.. code-block:: bash

    java -jar aacon-client.jar -i=d:\input.fasta -o=d:\result.txt -r="Quick conservation"

Get the list of presets supported by AACon web service

.. code-block:: bash

    java -jar aacon-client.jar -presets


------------

.. _web_cli:

Structure of the command line client
------------------------------------

+------------------------+--------------------------------------------------------------------------------------------+
| Packages               | Classes and Interfaces                                                                     |
+========================+============================================================================================+
| compbio.data.msa       | Annotation interface for all scoring methods                                               |
+------------------------+--------------------------------------------------------------------------------------------+
| compbio.data.sequence  | AAConWS data types                                                                         |
+------------------------+--------------------------------------------------------------------------------------------+
| compbio.metadata       | AAConWS meta data types                                                                    |
+------------------------+--------------------------------------------------------------------------------------------+
| compbio.ws.client      | AAConWS command line client package with AAConClient.java class containing the main method |
+------------------------+--------------------------------------------------------------------------------------------+

.. Please refer to the `data model javadoc`_ for a detailed description of each class and its methods.


------------

.. _web_functions:

AACon web service functions overview
------------------------------------


Functions for conservation calculation
::

    String id = analize(List<FastaSequence> list)
    String id = customAnalize(List<FastaSequence> sequenceList, List<Option> optionList)
    String id = presetAnalize(List<FastaSequence> sequenceList, Preset preset)

Functions pertaining to job monitoring and control
::

    JobStatus status = getJobStatus(String id)
    HashSet<Score> conservation = getConservation(String id)
    boolean cancelled = cancelJob(String id)
    ChunkHolder chunk = pullExecStatistics(String id, long marker)

Functions relating to service features discovery
::

    RunnerConfig rc = getRunnerOptions()
    Limit limit = getLimit(String name)
    LimitsManager lm = getLimits()
    PresetManager pm = getPresets()

.. Please refer to the `data model javadoc`_ for a detailed description of each class and its methods.

------------

.. _web_artifacts:

Building web services artifacts
-------------------------------


AAConWS are the standard `JAX-WS`_ SOAP web services, which are `WS-I`_ basic profile compatible. This means that you could use your favorite programming language to work with AAConWS. Below is how you can generate portable artifacts to work with AAConWS from Java. However, if programming in Java we recommend using our `client library`_ as it provides a handful of useful methods in addition to plain data types.

wsimport -keep http://www.compbio.dundee.ac.uk/aacon/AAConWS?wsdl


------------

.. _web_conn:

Connecting to AAConWS
---------------------


All the examples below assume that AACon web service command line client is in the classpath. You can download it from the download page. The code excerpt below will connect your program to AAConWS web service deployed in the University of Dundee.

::

    Annotation<AAConWS> client = AAConWSClient.connect();


If you want to work with the `generated artifacts`_ directly you can inspect AAConWS `command line client source code`_ and use it as a template for building your own custom AAConWS clients.


------------

.. _web_cons:

Calculating conservation
------------------------

Given that ``client`` is web service proxy, created as described in "Connecting to AAConWS" section, the conservation scores can be obtained as follows:

::

    1) List<FastaSequence> fsl = SequenceUtil.readFasta(new FileInputStream("alignment.fasta"));
    2) String jobId = client.analize(fsl);
    3) HashSet<Score> result = client.getAnnotation(jobId);


Line one loads sequence alignment from the file
Line two submits them to web service represented by AAConWS proxy
Line three retrieves the conservation scores calculated according to Shenkin algorithm. This line blocks the execution until the result is available. Use this with caution. In general, you should make sure that the calculation has been completed before attempting retrieving results. This is to avoid keeping the connection to the server on hold for a prolonged periods of time. While this may be ok with your local server, our public server (www.compbio.dundee.ac.uk/aacon) will not let you hold the connection for longer than 10 minutes. This is done to prevent excessive load on the server. The next section describes how to check the status of the calculation.
Methods and classes mentioned in the excerpt are available from the AAConWS client library.


------------

.. _web_status:

Checking the status of the calculation
--------------------------------------

You may have noticed that there was no pause between submitting the job and retrieving the results. This is because getAnnotation(jobId) method block the processing until the calculation is completed. However, taking into account that the connection holds server resources, our public server (www.compbio.dundee.ac.uk/aacon) is configured to reset the connection after 10 minutes of waiting. To work around the connection reset you are encouraged to check whether the calculation has been completed before accessing the results.	You can do it like this:

.. code-block:: java

    while (client.getJobStatus(jobId) != JobStatus.FINISHED) {
        Thread.sleep(2000); // wait two seconds, then recheck the status
    }


------------

.. _web_cons_presets:

Calculating conservation with presets
-------------------------------------

::

    1) PresetManager<AACon> presets = client.getPresets();
    2) String jobId = client.presetAnalize(fsl, presets.getPresetByName("Quick conservation"));
    3) HashSet<Score> result = client.getAnnotation(jobId);


Line one obtains the lists of presets supported by a web service.
Line two calls web service ``presetAnalise`` method with a chosen preset. This call returns a job identifier.
Lines three retrieves the results using job identifier.

Available presets are:

* "Quick conservation" (a collection of 12 fast conservation algorithms includes KABAT, JORES, SCHNEIDER, SHENKIN, GERSTEIN, TAYLOR_GAPS, TAYLOR_NO_GAPS, ZVELIBIL,ARMON, THOMPSON, NOT_LANCET, MIRNY, WILLIAMSON)
* "Slow conservation" (a collection of time consuming conservation algorithms includes LANDGRAF, KARLIN, SANDER, VALDAR and SMERFS)
* "Complete conservation" (all available algorithms)


------------

.. _web_cons_custom:

Calculating conservation with custom parameters
-----------------------------------------------

Below is the example of using custom parameters for SMERFS method.

::

    // Using options
    1) RunnerConfig<AACon> options = client.getRunnerOptions();
    2) options.getArgument("Calculation method").setDefaultValue("SMERFS");
    3) options.getArgument("SMERFS Column Scoring Method").setDefaultValue("MAX_SCORE");
    4) options.getArgument("SMERFS Gap Threshhold").setDefaultValue("1");
    5) String jobId = client.customAnalize(fsl, options.getArguments());
    6) HashSet<Score> result = client.getAnnotation(jobId);


Line one obtains the RunnerConfig object that holds information on supported parameters and their values
Line two retrieve a particular parameter from the holder by its name and sets the new value for this parameter.
Lines three and four do the same but for two more parameters
Line five submit a job to a web service
Line six retrieves the results of the analysis. The names of all the parameters supported by a web service can be obtained using options.getArguments() method. Further details on the methods available from RunnerConfig object are available from the javadoc.


------------

.. _web_cons_example:

A complete AAConWS web service client example
---------------------------------------------

Finally, a complete example of the program that connects to AAConWS service is below. The text notes are commented by block style comments e.g. /\* comment \*/, the code alternatives are commented out with the line comments - "//". You may want to remove line style comments to test alternatives of the functions. All you need for this to work is a AAConWS binary client. Please make sure that the client is in the Java class path before running this example.

.. code-block:: java

    import java.io.ByteArrayInputStream;
    import java.io.FileNotFoundException;
    import java.io.IOException;
    import java.util.List;
    import java.util.Set;

    import compbio.data.msa.Annotation;
    import compbio.data.sequence.FastaSequence;
    import compbio.data.sequence.Score;
    import compbio.data.sequence.SequenceUtil;
    import compbio.metadata.JobSubmissionException;
    import compbio.metadata.Preset;
    import compbio.metadata.PresetManager;
    import compbio.metadata.ResultNotAvailableException;
    import compbio.metadata.UnsupportedRuntimeException;
    import compbio.metadata.WrongParameterException;
    import compbio.runner.conservation.AACon;


    /**
     * AAConWS client example
     */
    public class AAConWSClientExample {

    	/*
    	 * Input sequences for alignment. For the simplicity keep them in the class
    	 */
    	static final String input = ">Foo      \r\n"
    	+ "MTADGPRELLQLRAAVRHRPQDFVAWLMLADAELGMGDTTAGEMAVQRGLALHPGHPEAV\r\n"
    	+ "ARLGRVRWTQQRHAEAAVLLQQASDAAPEHPGIALWLGHALEDAGQAEAAAAAYTRAHQL\r\n"
    	+ "LPEEPYITAQLLNWRRRLCDWRALDVLSAQVRAAVAQGVGAVEPFAFLSEDASAAEQLAC\r\n"
    	+ "ARTRAQAIAASVRPLAPTRVRSKGPLRVGFVSNGFGAHPTGLLTVALFEALQRRQPDLQM\r\n"
    	+ "HLFATSGDDGSTLRTRLAQASTLHDVTALGHLATAKHIRHHGIDLLFDLRGWGGGGRPEV\r\n"
    	+ "FALRPAPVQVNWLAYPGTSGAPWMDYVLGDAFALPPALEPFYSEHVLRLQGAFQPSDTSR\r\n"
    	+ "VVAEPPSRTQCGLPEQGVVLCCFNNSYKLNPQSMARMLAVLREVPDSVLWLLSGPGEADA\r\n"
    	+ "RLRAFAHAQGVDAQRLVFMPKLPHPQYLARYRHADLFLDTHPYNAHTTASDALWTGCPVL\r\n"
    	+ "TTPGETFAARVAGSLNHHLGLDEMNVADDAAFVAKAVALASDPAALTALHARVDVLRRES\r\n"
    	+ "GVFEMDGFADDFGALLQALARRHGWLGI\r\n"
    	+ "\r\n"
    	+ ">Bar                    \r\n"
    	+ "-----------------------------------MGDTTAGEMAVQRGLALH-------\r\n"
    	+ "---------QQRHAEAAVLLQQASDAAPEHPGIALWL-HALEDAGQAEAAAA-YTRAHQL\r\n"
    	+ "LPEEPYITAQLLN--------------------AVAQGVGAVEPFAFLSEDASAAE----\r\n"
    	+ "----------SVRPLAPTRVRSKGPLRVGFVSNGFGAHPTGLLTVALFEALQRRQPDLQM\r\n"
    	+ "HLFATSGDDGSTLRTRLAQASTLHDVTALGHLATAKHIRHHGIDLLFDLRGWGGGGRPEV\r\n"
    	+ "FALRPAPVQVNWLAYPGTSGAPWMDYVLGDAFALPPALEPFYSEHVLRLQGAFQPSDTSR\r\n"
    	+ "VVAEPPSRTQCGLPEQGVVLCCFNNSYKLNPQSMARMLAVLREVPDSVLWLLSGPGEADA\r\n"
    	+ "RLRAFAHAQGVDAQRLVFMPKLPHPQYLARYRHADLFLDTHPYNAHTTASDALWTGCPVL\r\n"
    	+ "TTPGETFAARVAGSLNHHLGLDEMNVADDAAFVAKAVALASDPAALTALHARVDVLRRES\r\n"
    	+ "GVFEMDGFADDFGALLQALARRHGWLGI\r\n"
    	+ "\r\n"
    	+ ">Noname             \r\n"
    	+ "-MTADGPRELLQLRAAVRHRPQDVAWLMLADAELGMGDTTAGEMAVQRGLALHPGHPEAV\r\n"
    	+ "ARLGRVRWTQQRHAEAAVLLQQASDAAPEHPGIALWLGHALED--------------HQL\r\n"
    	+ "LPEEPYITAQLDVLSAQVR-------------AAVAQGVGAVEPFAFLSEDASAAEQLAC\r\n"
    	+ "ARTRAQAIAASVRPLAPTRVRSKGPLRVGFVSNGFGAHPTGLLTVALFEALQRRQPDLQM\r\n"
    	+ "HLFATSGDDGSTLRTRLAQASTLHDVTALGHLATAKHIRHHGIDLLFDLRGWGGGGRPEV\r\n"
    	+ "FALRPAPVQVNWLAYPGTSGAPWMDYVLGDAFALPPALEPFYSEHVLRLQGAFQPSDTSR\r\n"
    	+ "VVAEPPSRTQCGLPEQGVVLCCFNNSYKLNPQSMARMLAVLREVPDSVLWLLSGPGEADA\r\n"
    	+ "RLRAFAHAQGVDAQRLVFMPKLPHPQYLARYRHADLFLDTHPYNAHTTASDALWTGCPVL\r\n"
    	+ "TTPGETFAARVAGSLNHHLGLDEMNVADDAAFVAKAVALASDPAALTALHARVDVLRRES\r\n"
    	+ "I---------------------------";

    	public static void main(String[] args) throws UnsupportedRuntimeException,
    			JobSubmissionException, WrongParameterException,
    			FileNotFoundException, IOException, ResultNotAvailableException,
    			InterruptedException {

    		/*
    		 * Annotation interface for AAConWS web service instance
    		 */
    		Annotation<AACon> client = AAConClient.connect();

    		/* Get the list of available presets */
    		PresetManager presetman = client.getPresets();

    		/* Get the Preset object by preset name */
    		Preset preset = presetman.getPresetByName("Complete conservation");

    		/*
    		 * Load sequences in FASTA format from the file You can use something
    		 * like new FileInputStream() to load sequence from the file
    		 */
    		List<FastaSequence> fastalist = SequenceUtil
    				.readFasta(new ByteArrayInputStream(input.getBytes()));

    		/*
    		 * Submit loaded sequences for an alignment using preset. The job
    		 * identifier is returned by this method, you can retrieve the results
    		 * with it sometime later.
    		 */
    		String jobId = client.presetAnalize(fastalist, preset);

    		/* This method will block for the duration of the calculation */
    		Set<Score> result = client.getAnnotation(jobId);

    		/*
    		 * This is a better way of obtaining results, it does not involve
    		 * holding the connection open for the duration of the calculation,
    		 * Besides, as the University of Dundee public server will reset the
    		 * connection after 10 minutes of idling, this is the only way to obtain
    		 * the results of long running task from our public server.
    		 */
    		// while (client.getJobStatus(jobId) != JobStatus.FINISHED) {
    		// Thread.sleep(1000); // wait a second, then recheck the status
    		// }

    		/* Output the alignment to standard out */
    		Score.write(result, System.out);

    		/* Alternatively, you can record retrieved alignment into the file */
    		// FileOutputStream out = new FileOutputStream("result.txt");
    		// Score.write(result, out);
    		// out.close();
    	}
    }


For a more detailed description of all available types and their functions please refer to the `data model javadoc`_.

.. links
.. _Java: http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html
.. _data model javadoc: ../docs/javadoc/index.html
.. _aacon homepage: ../index.html
.. _WS-I: http://www.ws-i.org/
.. _JAX-WS: http://jax-ws.java.net/
.. _client library: library.html
.. _generated artifacts: web_service.html#building_web_services_artifacts
.. _command line client source code: ../AAConClient.java
