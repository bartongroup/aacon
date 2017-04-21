Java Library
============

This page provides information pertaining to AACon Java library usage.


------------

.. _lib_usage:

Using the library
-----------------

.. note:: For impatient please find the complete example `below`_!

To calculate conservation using AACon follow the steps:

1. :ref:`lib_load` - Load data
2. :ref:`lib_init` - Initialise the executor
3. :ref:`lib_norm` - Decide on normalisation
4. :ref:`lib_meth` - Set conservation methods
5. :ref:`lib_output` - Output the results

All these steps are discussed in greater details below.


------------

.. _lib_load:

Loading data
------------

AACon methods require a List of FastaSequence objects as input. Below is an example of how to load the Fasta formatted sequence file or the Clustal alignment file. For this call to succeed all the sequences must be of the same length (as this is an alignment).

.. code-block:: java

    List<FastaSequence> sequences = CmdParser.openInputStream(<PATH_TO_INPUT_FILE>);


For the examples of the input file please see AACon `standalone executable help page`_.

------------

.. _lib_init:

Initialising the executor
-------------------------


AACon methods require the `ExecutorService`_ object which is used to parallel the calculations. The ExecutorService initialised with the number of threads equals to the number of cores on the executing machine is recommended for a faster calculation and can be obtained as follows:

.. code-block:: java

    int corenum = Runtime.getRuntime().availableProcessors();
    ExecutorService executor = Executors.newFixedThreadPool(corenum);

Please take care to initialise and pass only one executor to all the methods to avoid the waist of resources. After use, the executor must be disposed of, it can be done as follows:

.. code-block:: java

    executor.shutdown();


------------

.. _lib_norm:

Normalisation
-------------

AACon methods require the boolean parameter telling the system whether the results should be normalised or not. Normalised results have values between 0 and 1. Please note however, that some results cannot be normalised. In such a case, the system returns not normalised values, and log the issue to the standard error stream. The following formula is used for normalisation

.. code-block:: java

    n = (d - dmin)/(dmax - dmin)

Negative results first converted to positive by adding an absolute value of the most negative result.


------------

.. _lib_meth:

Setting the conservation methods
--------------------------------

AACon ``ConservationCalculator.getConservation()`` method takes a Set of ``ConservationMethod`` to use for calculating the conservation. Below is a few examples of making such sets.

* To calculate conservation according to KABAT method construct a set in the following way -

    .. code-block:: java

        EnumSet.of(ConservationMethod.KABAT).


* For all the methods use the following construct -

    .. code-block:: java

        EnumSet.allOf(ConservationMethod.class)


* For a set of methods including KABAT, JORES, SCHNEIDER, SHENKIN and GERSTEIN use the following construct -

    .. code-block:: java

        EnumSet.range(ConservationMethod.KABAT, ConservationMethod.GERSTEIN)


------------

.. _lib_output:

Output the results
------------------


The results can be output using one of the following methods:

1. *ConservationFormatter.formatResults(Map<ConservationMethod, double[]> scores, OutputStream outStream)*
2. *ConservationFormatter.formatResults(Map<ConservationMethod, double[]> scores, String outFilePath, Format format, List<FastaSequence> alignment)*

Use the first method to output the results of the calculation without an alignment to any OutputStream. Use the second method to output results with the alignment.

First method usage example:

.. code-block:: java

    // Usage example - printing results to the console
    ConservationFormatter.outputScoreLine(result, System.out)
    // printing results to the file
    FileOutputStream outfile = new FileOutputStream("results.txt");
    ConservationFormatter.formatResults(result, outfile);
    outfile.close();

Second method usage example:

.. code-block:: java

    // Usage example - printing results with alignment in Fasta format to the file called output.txt
    ConservationFormatter.formatResults(result, "test.txt", Format.RESULT_WITH_ALIGNMENT, sequences);


For more information on the output formats please see standalone `AACon help`_.


------------

.. _lib_example:

Example
-------

Using AACon library for calculating conservation:

.. code-block:: java

    // Determine the number of CPU cores available on the system.
    int corenum = Runtime.getRuntime().availableProcessors();
    // Initialize the Executor instance with a number of cores
    ExecutorService executor = Executors.newFixedThreadPool(corenum);

    // Load the data from the file containing either Clustal formatted alignment
    // or a list of FASTA formatted sequences. Assuming that small.align file is
    // in the same directory as this program
    List<FastaSequence> sequences = CmdParser.openInputStream("small.align");

    // Calculate conservation scores using all methods.
    Map<ConservationMethod, double[]> result = getConservation(sequences, true,
                             EnumSet.allOf(ConservationMethod.class)), executor);

    // Print the results to the console.
    ConservationFormatter.formatResults(result, System.out);


.. links
.. _below: library.html#example
.. _standalone executable help page: client.html
.. _ExecutorService: http://download.oracle.com/javase/1.5.0/docs/api/java/util/concurrent/ExecutorService.html
.. _AACon help: client.html]#the_output
