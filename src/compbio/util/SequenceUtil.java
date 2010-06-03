/* Copyright (c) 2009 Peter Troshin
 *  
 *  Jalview Web Services @version: 2.0     
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

package compbio.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Utility class for operations on sequences
 * 
 * @author pvtroshin
 * 
 *         Date September 2009
 */
public final class SequenceUtil {

    /**
     * A whitespace character: [\t\n\x0B\f\r]
     */
    public static final Pattern WHITE_SPACE = Pattern.compile("\\s");

    /**
     * A digit
     */
    public static final Pattern DIGIT = Pattern.compile("\\d");

    /**
     * Non word
     */
    public static final Pattern NONWORD = Pattern.compile("\\W");

    /**
     * Valid Amino acids
     */
    public static final Pattern AA = Pattern.compile("[ARNDCQEGHILKMFPSTWYV]+",
	    Pattern.CASE_INSENSITIVE);

    /**
     * inversion of AA pattern
     */
    public static final Pattern NON_AA = Pattern.compile(
	    "[^ARNDCQEGHILKMFPSTWYV]+", Pattern.CASE_INSENSITIVE);

    /**
     * Same as AA pattern but with two additional letters - XU
     */
    public static final Pattern AMBIGUOUS_AA = Pattern.compile(
	    "[ARNDCQEGHILKMFPSTWYVXU]+", Pattern.CASE_INSENSITIVE);

    /**
     * Nucleotides a, t, g, c, u
     */
    public static final Pattern NUCLEOTIDE = Pattern.compile("[AGTCU]+",
	    Pattern.CASE_INSENSITIVE);

    /**
     * Ambiguous nucleotide
     */
    public static final Pattern AMBIGUOUS_NUCLEOTIDE = Pattern.compile(
	    "[AGTCRYMKSWHBVDNU]+", Pattern.CASE_INSENSITIVE); // see IUPAC
    /**
     * Non nucleotide
     */
    public static final Pattern NON_NUCLEOTIDE = Pattern.compile("[^AGTCU]+",
	    Pattern.CASE_INSENSITIVE);

    private SequenceUtil() {
    } // utility class, no instantiation

    /*
     * public static void write_PirSeq(OutputStream os, FastaSequence seq)
     * throws IOException { BufferedWriter pir_out = new BufferedWriter(new
     * OutputStreamWriter(os)); pir_out.write(">P1;" + seq.getId() +
     * SysPrefs.newlinechar); pir_out.write(seq.getSequence() +
     * SysPrefs.newlinechar); pir_out.close(); }
     * 
     * public static void write_FastaSeq(OutputStream os, FastaSequence seq)
     * throws IOException { BufferedWriter fasta_out = new BufferedWriter( new
     * OutputStreamWriter(os)); fasta_out.write(">" + seq.getId() +
     * SysPrefs.newlinechar); fasta_out.write(seq.getSequence() +
     * SysPrefs.newlinechar); fasta_out.close(); }
     */

    /**
     * @return true is the sequence contains only letters a,c, t, g, u
     */
    public static boolean isNucleotideSequence(FastaSequence s) {
	return isNonAmbNucleotideSequence(s.getSequence());
    }

    /**
     * Ambiguous DNA chars : AGTCRYMKSWHBVDN // differs from protein in only one
     * (!) - B char
     */
    public static boolean isNonAmbNucleotideSequence(String sequence) {
	sequence = cleanSequence(sequence);
	if (DIGIT.matcher(sequence).find()) {
	    return false;
	}
	if (NON_NUCLEOTIDE.matcher(sequence).find()) {
	    return false;
	    /*
	     * System.out.format("I found the text starting at " +
	     * "index %d and ending at index %d.%n", nonDNAmatcher .start(),
	     * nonDNAmatcher.end());
	     */
	}
	Matcher DNAmatcher = NUCLEOTIDE.matcher(sequence);
	return DNAmatcher.find();
    }

    /**
     * Removes all whitespace chars in the sequence string
     * 
     * @param sequence
     * @return cleaned up sequence
     */
    public static String cleanSequence(String sequence) {
	assert sequence != null;
	final Matcher m = WHITE_SPACE.matcher(sequence);
	sequence = m.replaceAll("").toUpperCase();
	return sequence;
    }

    /**
     * Removes all special characters and digits as well as whitespace chars
     * from the sequence
     * 
     * @param sequence
     * @return cleaned up sequence
     */
    public static String deepCleanSequence(String sequence) {
	sequence = cleanSequence(sequence);
	sequence = DIGIT.matcher(sequence).replaceAll("");
	sequence = NONWORD.matcher(sequence).replaceAll("");
	Pattern othernonSeqChars = Pattern.compile("[_-]+");
	sequence = othernonSeqChars.matcher(sequence).replaceAll("");
	return sequence;
    }

    /**
     * 
     * @param sequence
     * @return true is the sequence is a protein sequence, false overwise
     */
    public static boolean isProteinSequence(String sequence) {
	sequence = cleanSequence(sequence);
	if (isNonAmbNucleotideSequence(sequence)) {
	    return false;
	}
	if (DIGIT.matcher(sequence).find()) {
	    return false;
	}
	if (NON_AA.matcher(sequence).find()) {
	    return false;
	}
	Matcher protmatcher = AA.matcher(sequence);
	return protmatcher.find();
    }

    /**
     * Check whether the sequence confirms to amboguous protein sequence
     * 
     * @param sequence
     * @return return true only if the sequence if ambiguous protein sequence
     *         Return false otherwise. e.g. if the sequence is non-ambiguous
     *         protein or DNA
     */
    public static boolean isAmbiguosProtein(String sequence) {
	sequence = cleanSequence(sequence);
	if (isNonAmbNucleotideSequence(sequence)) {
	    return false;
	}
	if (DIGIT.matcher(sequence).find()) {
	    return false;
	}
	if (NON_AA.matcher(sequence).find()) {
	    return false;
	}
	if (AA.matcher(sequence).find()) {
	    return false;
	}
	Matcher amb_prot = AMBIGUOUS_AA.matcher(sequence);
	return amb_prot.find();
    }

    /**
     * Writes list of FastaSequeces into the outstream formatting the sequence
     * so that it contains width chars on each line
     * 
     * @param outstream
     * @param sequences
     * @param width
     *            - the maximum number of characters to write in one line
     * @throws IOException
     */
    public static void writeFasta(OutputStream outstream,
	    List<FastaSequence> sequences, int width) throws IOException {
	OutputStreamWriter writer = new OutputStreamWriter(outstream);
	BufferedWriter fastawriter = new BufferedWriter(writer);
	for (FastaSequence fs : sequences) {
	    fastawriter.write(fs.getOnelineFasta());
	}
	outstream.flush();
	fastawriter.close();
	writer.close();
    }

    /**
     * Reads fasta sequences from inStream into the list of FastaSequence
     * objects
     * 
     * @param inStream
     *            from
     * @return list of FastaSequence objects
     * @throws IOException
     */
    public static List<FastaSequence> readFasta(InputStream inStream)
	    throws IOException {
	List<FastaSequence> seqs = new ArrayList<FastaSequence>();

	BufferedReader infasta = new BufferedReader(new InputStreamReader(
		inStream, "UTF8"), 16000);
	Pattern pattern = Pattern.compile("//s+");

	String line;
	String sname = "", seqstr = null;
	do {
	    line = infasta.readLine();
	    if (line == null || line.startsWith(">")) {
		if (seqstr != null)
		    seqs.add(new FastaSequence(sname.substring(1), seqstr));
		sname = line; // remove >
		seqstr = "";
	    } else {
		String subseq = pattern.matcher(line).replaceAll("");
		seqstr += subseq;
	    }
	} while (line != null);

	infasta.close();
	return seqs;
    }

    /**
     * Writes FastaSequence in the file, each sequence will take one line only
     * 
     * @param os
     * @param sequences
     * @throws IOException
     */
    public static void writeFasta(OutputStream os, List<FastaSequence> sequences)
	    throws IOException {
	OutputStreamWriter outWriter = new OutputStreamWriter(os);
	BufferedWriter fasta_out = new BufferedWriter(outWriter);
	for (FastaSequence fs : sequences) {
	    fasta_out.write(fs.getOnelineFasta());
	}
	fasta_out.close();
	outWriter.close();
    }

}
