Standalone Client
=================

This page provides some tips on using AACon command line tool.

------------

.. _cli_exec:

Executing AACon
---------------

Type java -jar <AACon jar file name>

.. code-block:: bash

    java -jar compbio-conservation-1.1.jar

This will print AACon help.


------------

.. _cli_input:

Input
-----

Input to AACon can be either

* The list of aligned FASTA formatted sequences
* The Clustal formatted alignment

**FASTA example**

::

    >UniRef90_Q5VFX
    MEWLKIALVAFSGRGNDVDYPCSILGRVANNDKELVNILRNGFFLLTYRMNFSPLPHSSVTSDKGWGCLVRS
    SQMLLAHALWRYSANDCRLDHFRDMDTEDSTPFSLHKMVRAVMKKADVFRPEYWTPSQGCEAIRCCVNNAVD
    RKLIPPIRVVVCSQGCLLAREICSNLEFGTVLILAPMRCGASRRMTQMMFFSLEHLLHSSACIGVVGGVPQR
    SYYILGTSGQRLLYLDPHCMTQEALVSSHAEKAGVVTVTASLVKSVRWDCVDTSCFLGFLVDSFAEWLELRT
    CLEELQRRGMEQLLCVDDGVAVGLVDEEIAGWPSEEDVAE
    >UniRef90_Q4VBK1
    ---------------PDTDEPVWILGACYNVKTKKSELLSDSRLWFTYRKKFSPIGGTGPSSDAGWGCMLRC
    GQMILAQALWRWDPEKHQPKEYQRFLDKKDSCYSIHQMAQMGVGE-GKSVGEWYGPNTVAQVLKKL----AL
    FDDWNSLSVYVSMDNTVVIEDIKKLCDWRPLLLVIPLRMGINS-INPVYIQALKECFKMPQSCGVLGGKPNL
    AYYFIGFIDDELIYLDPHTTQQA-VDTESGSAVDDQSFHCQRPHRMKITSLDPSVALGFFCKSEEDFDSWCD
    LVQ-------QELLKKRNLRMFELVEKHPSHWPPF-----
    >UniRef90_A7RTH8
    ---MDAACVTYEGVSHETEEDVWILGKRYNIDMGYLNTDVRSRIWLTYRKNFPKIGGTGPTTDSGWGCMLRC
    GQMMLAQALWQWDPENEYMQILEAFLDKKDSLYSIHQIAQMGVSE-GKAVGSWFGPNTVAQVLKKL----SA
    FDDWSSLCLHVAMDNTVIIEDISN---WRPLVLFIPLRLGLTEM-NVVYNEPLKACFTFKQSLGIIGGRPNH
    ATYFIGYFGNNLVYLDPHTTQQTVN-PDELSRIPDGSFHCVYPCRMNIADVDPSVALGFFCKSEEDFDDLCQ
    QIKKIIDGKSRPMFEIAK--------DRPQHWPVLE----


**Clustal example**

::

    CLUSTAL

    QUERY/1-328 MEWLKIALVAFSGRGNDVDYPCSILGRVANNDKELVNILRNGFFLLTYRMNFSPLPHSSV
    UniRef90_Q4VBK1/1-294 ---------------PDTDEPVWILGACYNVKTKKSELLSDSRLWFTYRKKFSPIGGTGP
    UniRef90_A7RTH8/1-303 ---MDAACVTYEGVSHETEEDVWILGKRYNIDMGYLNTDVRSRIWLTYRKNFPKIGGTGP

    QUERY/1-328 TSDKGWGCLVRSSQMLLAHALWRYSANDCRLDHFRDMDTEDSTPFSLHKMVRAVMKKADV
    UniRef90_Q4VBK1/1-294 SSDAGWGCMLRCGQMILAQALWRWDPEKHQPKEYQRFLDKKDSCYSIHQMAQMGVGEGKS
    UniRef90_A7RTH8/1-303 TTDSGWGCMLRCGQMMLAQALWQWDPENEYMQILEAFLDKKDSLYSIHQIAQMGVSEGKA

    QUERY/1-328 FRPEYWTPSQGCEAIRCCVNNAVDRKLIPPIRVVVCSQGCLLAREICSNLEFGTVLILAP
    UniRef90_Q4VBK1/1-294 VG-EWYGPN----TVAQVLKKLALFDDWNSLSVYVSMDNTVVIEDIKKLCDWRPLLLVIP
    UniRef90_A7RTH8/1-303 VG-SWFGPN----TVAQVLKKLSAFDDWSSLCLHVAMDNTVIIEDIS---NWRPLVLFIP

    QUERY/1-328 MRCGASRRMTQMMFFSLEHLLHSSACIGVVGGVPQRSYYILGTSGQRLLYLDPHCMTQEA
    UniRef90_Q4VBK1/1-294 LRMGIN-SINPVYIQALKECFKMPQSCGVLGGKPNLAYYFIGFIDDELIYLDPH-TTQQA
    UniRef90_A7RTH8/1-303 LRLGLT-EMNVVYNEPLKACFTFKQSLGIIGGRPNHATYFIGYFGNNLVYLDPH-TTQQT

    QUERY/1-328 LVSSHAEKAGVVTVTASLVKSVRWDCVDTSCFLGFLVDSFAEWLELRTCLEELQRRGMEQ
    UniRef90_Q4VBK1/1-294 VDTESGSAVDDQSFHCQRPHRMKITSLDPSVALGFFCKSEEDFDSWCDLVQQ-------E
    UniRef90_A7RTH8/1-303 VNPDELSRIPDGSFHCVYPCRMNIADVDPSVALGFFCKSEEDFDDLCQQIK--------K

    QUERY/1-328 LLCVDDGVAVGLVDEEIAGWPSEEDVAE
    UniRef90_Q4VBK1/1-294 LLKKRNLRMFELVEKHPSHWPPF-----
    UniRef90_A7RTH8/1-303 IIDGKSRPMFEIAKDRPQHWPVLE----


The gaps in the alignment can be represented by \*, -, space character, X and . (a dot). Other, custom gap characters can also be defined.


------------

.. _cli_output:

Output
------

AACon supports two modes of the output the first and the default mode is *RESULT_NO_ALIGNMENT*. Only the method names and the conservation score for each position in the alignment are presented. The second output format is *RESULT_WITH_ALIGNMEN*, where the input, in the form of FASTA formatted alignment, is included at the beginning of the file, followed by the same output as *RESULT_NO_ALIGNMENT*.

RESULT_NO_ALIGNMENT example is below:

::

    #KABAT 0.75 0.75 0.75 0.375 0.375 0.375 0.938 0.375 0.938 0.375 0.37 0.3
    #JORES 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    #SCHNEIDER 0.667 0.667 0.667 0.333 0.333 0.333 0.754 0.333 0.754 0.333 0
    #SHENKIN 0.779 0.779 0.779 0.46 0.46 0.46 0.845 0.46 0.845 0.46 0.46 0.4
    #GERSTEIN 0.667 0.667 0.667 0.333 0.333 0.333 0.754 0.333 0.754 0.33 0.3
    #TAYLOR_GAPS 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.632 1 0 1 0.526 0.842 0 0.5
    #TAYLOR_NO_GAPS 0 1 0.526 0 0.895 0 0.947 0 0 0.789 0 0.842 0.947 0 0.94
    #ZVELIBIL 0.1 0.3 0.4 0.1 0.2 0.1 0.3 0.1 0.3 0.2 0.2 0.1 0.3 0 0.2 0.2
    #KARLIN 0.338 0.338 0.338 0.433 0.31 0.264 0.558 0.288 0.558 0.338 0.431
    #ARMON 0.617 0.617 0.617 0.208 0.103 0.063 0.617 0.129 0.617 0.177 0.204
    #THOMPSON 0.719 0.727 0.708 0.451 0.394 0.384 0.471 0.374 0.49 0.428 0.4
    #NOT_LANCET 0.011 0.005 0.054 0.172 0.14 0.108 0.194 0.145 0.194 0.1 0.2
    #MIRNY 0.667 0.667 0.667 0.754 0.333 0.754 0.754 0.754 0.754 0.333 0.754
    #WILLIAMSON 0.884 0.953 0.971 0.503 0.511 0.445 0.645 0.552 0.503 0.522
    #LANDGRAF 0.752 0.726 0.847 0.539 0.436 0.25 0.433 0.476 0.547 0.308 0.6
    #SANDER 0 0 0 0.16 0.107 0.093 0.133 0.027 0.16 0.12 0.2 0.107 0.17 0.01
    #VALDAR 0 0 0 0.187 0.115 0.115 0.256 0.046 0.256 0.164 0.233 0.092 0.25
    #SMERFS 0.625 0.625 0.625 0.625 0.446 0.3 0.232 0.232 0.232 0.232 0.232


RESULT_WITH_ALIGNMENT output example is below:

::

    >QUERY
    MEWLKIALVAFSGRGNDVDYPCSILGRVANNDKELVNILRNGFFLLTYRMNFSPLPHSSVTSDKGWGCLVRSSQMLLAHA
    LWRYSANDCRLDHFRDMDTEDSTPFSLHKMVRAVMKKADVFRPEYWTPSQGCEAIRCCVNNAVDRKLIPPIRVVVCSQGC
    LLAREICSNLEFGTVLILAPMRCGASRRMTQMMFFSLEHLLHSSACIGVVGGVPQRSYYILGTSGQRLLYLDPHCMTQEA
    LVSSHAEKAGVVTVTASLVKSVRWDCVDTSCFLGFLVDSFAEWLELRTCLEELQRRGMEQLLCVDDGVAVGLVDEEIAGW
    PSEEDVAE
    >UniRef90_Q4VBK1
    ---------------PDTDEPVWILGACYNVKTKKSELLSDSRLWFTYRKKFSPIGGTGPSSDAGWGCMLRCGQMILAQA
    LWRWDPEKHQPKEYQRFLDKKDSCYSIHQMAQMGVGE-GKSVGEWYGPNTVAQVLKKL----ALFDDWNSLSVYVSMDNT
    VVIEDIKKLCDWRPLLLVIPLRMGINS-INPVYIQALKECFKMPQSCGVLGGKPNLAYYFIGFIDDELIYLDPHTTQQA-
    VDTESGSAVDDQSFHCQRPHRMKITSLDPSVALGFFCKSEEDFDSWCDLVQ-------QELLKKRNLRMFELVEKHPSHW
    PPF-----
    >UniRef90_A7RTH8
    ---MDAACVTYEGVSHETEEDVWILGKRYNIDMGYLNTDVRSRIWLTYRKNFPKIGGTGPTTDSGWGCMLRCGQMMLAQA
    LWQWDPENEYMQILEAFLDKKDSLYSIHQIAQMGVSE-GKAVGSWFGPNTVAQVLKKL----SAFDDWSSLCLHVAMDNT
    VIIEDISN---WRPLVLFIPLRLGLTEM-NVVYNEPLKACFTFKQSLGIIGGRPNHATYFIGYFGNNLVYLDPHTTQQTV
    N-PDELSRIPDGSFHCVYPCRMNIADVDPSVALGFFCKSEEDFDDLCQQIKKIIDGKSRPMFEIAK--------DRPQHW
    PVLE----
    #KABAT 3 3 3 6 6 6 1.5 6 1.5 6 6 6 1.5 6 6 9 3 3 3 3 3 3 3 1 1 1 9 9 3 1 9 3 9 9 9 9 3 9 3 9 9 3 3 6 6 6 6 6 9 9 3 3 9 9 9 9 6 6 6 6 6 1.5 1.5 6 9 9 3 9 3 1 1 9 9 1.5 3 3 3 3
    #SANDER -72 -72 -72 -36 -48 -51 -42 -66 -36 -45 -27 -48 -33 -54 -45 3 30 9 30 -12 12 0 39 45 54 45 -18 -36 -72 -72 -72 -72


.. warning:: Please note that only a small part of the output is presented in all examples that follow for the sake of simplicity. Therefore the results are for illustration only. Normally each position of the alignment is given a conservation score by each calculation method.


``-f`` option can be used to specify the output format. For example the command

.. code-block:: bash

    java -jar compbio-conservation-1.1.jar -i=data.align -m=KABAT -f=RESULT_WITH_ALIGNMENT

outputs the input alignment into the result file.


------------

.. _cli_calcon:

Calculating conservation
------------------------


**Calculate conservation using KABAT conservation method**

The alignment is read from the data.align file.

.. code-block:: bash

    java -jar compbio-conservation-1.1.jar -i=data.align -m=KABAT


The above command will print the results to the console. You should see something like this:

::

    #KABAT 3 3 3 6 6 6 1.5 6 1.5 6 6 6 1.5 6 6 9 3 3 3 3 3 3 3 1 1 1 9 9


**Calculate conservation using two or more methods, for example SANDER and KABAT**

Just specify the list of comma separated methods after ``-m`` switch like this:

java -jar compbio-conservation-1.1.jar -i=data.align -m=SANDER,KABAT

The results will look something like this:
::

    #KABAT 3 3 3 6 6 6 1.5 6 1.5 6 6 6 1.5 6 6 9 3 3 3 3 3 3 3 1 1 1 9 9 3
    #SANDER -72 -72 -72 -36 -48 -51 -42 -66 -36 -45 -27 -48 -33 -54 -45 3 30


**Calculate conservation using all supported conservation methods and make results comparable**

.. code-block:: bash

    java -jar compbio-conservation-1.1.jar -i=data.align -n

The results will be printed to the console. Where ``-n`` normalizes all results.


------------

.. _cli_custom_gap:

Custom gap character
--------------------


Assuming that the gaps in the alignment are represented by *'-'*, *'_'* and *'x'* symbols the following command will interpret them correctly.

.. code-block:: bash

    java -jar compbio-conservation-1.1.jar -i=data.align -m=KABAT -g=-,_,x


------------

.. _cli_smerfs_custom:

Running SMERFS with custom parameters
-------------------------------------

Unlike other methods, SMERFS supports a few custom parameters, in particular

1. The window width - an integer and an odd number
2. Two methods of window scores to columns allocation:

  *MID_SCORE* - gives the window score to the middle column
  *MAX_SCORE* - gives the column the highest score of all the windows it belongs to

3. A gap percentage cutoff - a float greater than 0 and smaller or equal 1

For example:

.. code-block:: bash

    java -jar compbio-conservation-1.1.jar -i=data.align -m=SMERFS -s=5,MID_SCORE,0.1


------------

.. _cli_file_output:

Outputing to a file
-------------------

Use -o option to print the results to the file instead of a console. For example

.. code-block:: bash

    java -jar compbio-conservation-1.1.jar -i=data.align -n -o=outfile.txt


will produce output.txt results file. Nothing will be printed to the console.


------------

.. _cli_exec_details:

Execution details
-----------------

Use -d option to tell AACon to output its execution details.

.. code-block:: bash

    java -jar compbio-conservation-1.1.jar -i=data.align -n -d=stat.out -f=RESULT_WITH_ALIGNMENT

The output of the previous command, the context of the stat.out file is below

::

    No methods are request assuming all are required.
    No output file is provided, writing results to the standard output.
    Setting output format to RESULT_WITH_ALIGNMENT
    Using 4 CPUs
    Start time: 2010/12/08 12:02:35
    Alignment loaded in: 115 ms
    Alignment has: 3 sequences.
    Alignment length is: 328
    KARLIN 15 ms
    LANDGRAF 25 ms
    SANDER 2 ms
    VALDAR 2 ms
    SMERFS 3 ms
    SCHNEIDER 1 ms
    KABAT 1 ms
    SHENKIN 3 ms
    JORES 5 ms
    GERSTEIN 3 ms
    ARMON 1 ms
    THOMPSON 4 ms
    NOT_LANCET 4 ms
    ZVELIBIL 13 ms
    MIRNY 4 ms
    TAYLOR_GAPS 20 ms
    WILLIAMSON 4 ms
    TAYLOR_NO_GAPS 21 ms
    Total calculation time: 0 s
    End time: 2010/12/08 12:02:36


------------

.. _cli_results_norm:

Results normalization
---------------------

Different conservation algorithms produce vastly different values. To compare them meaningfully one need to bring the results into the same range. AACon bring the results to the range between 0 and 1 if ``-n`` switch is used.


------------

.. _cli_large_aln:

Conservation for large alignments
---------------------------------

.. attention:: Troubleshooting Java VM running out of memory

To enable AACon to deal with large alignments (thousands of sequences) let the Java Virtual Machine use more memory with the flag ``-Xmx<AMOUNT_OF_MEMORY>`` to your command line as follows:

.. code-block:: bash

    java -Xmx1G -jar aacon.jar <options>

Where ``1G = 1`` gigabyte of memory, the same result can be achieved with the following instruction ``-Xmx1000M``. However, floating point numbers like ``-Xmx1.5G`` are not allowed.
