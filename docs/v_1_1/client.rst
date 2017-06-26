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

RESULT_NO_ALIGNMENT (and normalised) example is below:

::

  #KABAT          0.7500 0.7500 0.7500 0.3750 0.3750 0.3750 0.9375 0.3750 0.9375 0.3750 0.3750 0.3750  ...
  #JORES          0.0000 0.0000 0.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000  ...
  #SCHNEIDER      0.6668 0.6668 0.6668 0.3332 0.3332 0.3332 0.7540 0.3332 0.7540 0.3332 0.3332 0.3332  ...
  #SHENKIN        0.7789 0.7789 0.7789 0.4600 0.4600 0.4600 0.8448 0.4600 0.8448 0.4600 0.4600 0.4600  ...
  #GERSTEIN       0.6667 0.6667 0.6667 0.3333 0.3333 0.3333 0.7540 0.3333 0.7540 0.3333 0.3333 0.3333  ...
  #TAYLOR_GAPS    0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000  ...
  #TAYLOR_NO_GAPS 0.0000 1.0000 0.5263 0.0000 0.8947 0.0000 0.9474 0.0000 0.0000 0.7895 0.0000 0.8421  ...
  #ZVELIBIL       0.1000 0.3000 0.4000 0.1000 0.2000 0.1000 0.3000 0.1000 0.3000 0.2000 0.2000 0.1000  ...
  #KARLIN         0.3377 0.3377 0.3377 0.4333 0.3099 0.2641 0.5584 0.2877 0.5584 0.3377 0.4308 0.3377  ...
  #ARMON          0.6171 0.6171 0.6171 0.2080 0.1034 0.0625 0.6171 0.1289 0.6171 0.1768 0.2036 0.1027  ...
  #THOMPSON       0.7195 0.7269 0.7076 0.4506 0.3942 0.3841 0.4710 0.3744 0.4896 0.4283 0.4041 0.3928  ...
  #NOT_LANCET     0.0108 0.0054 0.0538 0.1720 0.1398 0.1075 0.1935 0.1452 0.1935 0.1344 0.2097 0.1344  ...
  #MIRNY          0.6667 0.6667 0.6667 0.7540 0.3333 0.7540 0.7540 0.7540 0.7540 0.3333 0.7540 0.3333  ...
  #WILLIAMSON     0.8844 0.9531 0.9712 0.5029 0.5109 0.4447 0.6453 0.5525 0.5029 0.5220 0.6764 0.5195  ...
  #LANDGRAF       0.7531 0.7268 0.8473 0.5329 0.4288 0.2384 0.4247 0.4698 0.5405 0.2981 0.6343 0.3350  ...
  #SANDER         0.0000 0.0000 0.0000 0.1600 0.1067 0.0933 0.1333 0.0267 0.1600 0.1200 0.2000 0.1067  ...
  #VALDAR         0.0000 0.0000 0.0000 0.1867 0.1145 0.1145 0.2530 0.0482 0.2530 0.1627 0.2349 0.0904  ...
  #SMERFS         0.6250 0.6250 0.6250 0.6250 0.4462 0.3005 0.2324 0.2324 0.2324 0.2324 0.2324 0.2324  ...


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
    #KABAT          3.0000 3.0000 3.0000 6.0000 6.0000 6.0000 1.5000 6.0000 1.5000 6.0000 6.0000 6.0000 ...
    (...)
    #SMERFS         0.2593 0.2593 0.2593 0.2593 -0.0938 -0.3816 -0.5161 -0.5161 -0.5161 -0.5161 -0.5161 ...


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

    #KABAT          3.0000 3.0000 3.0000 6.0000 6.0000 6.0000 1.5000 6.0000 1.5000 6.0000 6.0000 6.0000 ...


**Calculate conservation using two or more methods, for example SANDER and KABAT**

Just specify the list of comma separated methods after ``-m`` switch like this:

java -jar compbio-conservation-1.1.jar -i=data.align -m=SANDER,KABAT

The results will look something like this:
::

    #KABAT          3.0000 3.0000 3.0000 6.0000 6.0000 6.0000 1.5000 6.0000 1.5000 6.0000 6.0000 6.0000 ...
    #SANDER         -72.0000 -72.0000 -72.0000 -36.0000 -48.0000 -51.0000 -42.0000 -66.0000 -36.0000 -45.0000 -27.0000 -48.0000 ...


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
    Using 8 CPUs
    Start time: 2017/05/23 09:58:13
    Alignment loaded in: 68 ms
    Alignment has: 3 sequences.
    Alignment length is: 328
    KARLIN 17 ms
    LANDGRAF 17 ms
    SANDER 2 ms
    VALDAR 4 ms
    SMERFS 4 ms
    KABAT 3 ms
    SHENKIN 3 ms
    SCHNEIDER 5 ms
    ARMON 3 ms
    JORES 6 ms
    GERSTEIN 8 ms
    THOMPSON 6 ms
    NOT_LANCET 4 ms
    MIRNY 9 ms
    WILLIAMSON 11 ms
    TAYLOR_NO_GAPS 20 ms
    ZVELIBIL 21 ms
    TAYLOR_GAPS 23 ms
    Total calculation time: 0 s
    End time: 2017/05/23 09:58:13


------------

.. _cli_results_norm:

Results normalization
---------------------

Different conservation algorithms produce vastly different values. To compare them meaningfully one need to bring the results into the same range. AACon bring the results to the range between 0 and 1 if ``-n`` switch is used.

.. warning:: Normalization should be used carefully, especially when working with very 'gappy' alignments, where the occupancy of a  single amino acid is observed in aligned columns.


------------

.. _cli_large_aln:

Conservation for large alignments
---------------------------------

.. .. attention:: Troubleshooting Java VM running out of memory

To enable AACon to deal with large alignments (thousands of sequences) let the Java Virtual Machine use more memory with the flag ``-Xmx<AMOUNT_OF_MEMORY>`` to your command line as follows:

.. code-block:: bash

    java -Xmx1G -jar aacon.jar <options>

Where ``1G = 1`` gigabyte of memory, the same result can be achieved with the following instruction ``-Xmx1000M``. However, floating point numbers like ``-Xmx1.5G`` are not allowed.
