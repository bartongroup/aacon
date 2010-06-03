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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;

/**
 * A FASTA formatted sequence. Please note that this class does not make any
 * assumptions as to what sequence it store e.g. it could be nucleotide, protein
 * or even gapped alignment sequence! The only guarantee it makes is that the
 * sequence does not contain white space characters e.g. spaces, new lines etc
 * 
 * @author pvtroshin
 * 
 *         Date September 2009
 */

@XmlAccessorType(XmlAccessType.FIELD)
public class FastaSequence {

    /**
     * Sequence id
     */
    private String id;

    // TODO what about gapped sequence here! should be indicated
    /**
     * Returns the string representation of sequence
     */
    private String sequence;

    private FastaSequence() {
	// Default constructor for JaxB
    }

    /**
     * Upon construction the any whitespace characters are removed from the
     * sequence
     * 
     * @param id
     * @param sequence
     */
    public FastaSequence(String id, String sequence) {
	this.id = id.trim();
	this.sequence = SequenceUtil.cleanSequence(sequence);
    }

    /**
     * Gets the value of id
     * 
     * @return the value of id
     */
    public String getId() {
	return this.id;
    }

    /**
     * Gets the value of sequence
     * 
     * @return the value of sequence
     */
    public String getSequence() {
	return this.sequence;
    }

    public static int countMatchesInSequence(final String theString,
	    final String theRegExp) {
	final Pattern p = Pattern.compile(theRegExp);
	final Matcher m = p.matcher(theString);
	int cnt = 0;
	while (m.find()) {
	    cnt++;
	}
	return cnt;
    }

    public String getFormattedFasta() {
	return getFormatedSequence(80);
    }

    /**
     * 
     * @return one		// TODO Auto-generated constructor stub line name, next line sequence, no matter what the sequence
     *         length is
     */
    public String getOnelineFasta() {
	String fasta = ">" + getId() + "\n";
	fasta += getSequence() + "\n";
	return fasta;
    }

    /**
     * Format sequence per width letter in one string. Without spaces.
     * 
     * @return multiple line formated sequence, one line width letters length
     * 
     */
    public String getFormatedSequence(final int width) {
	if (sequence == null) {
	    return "";
	}

	assert width >= 0 : "Wrong width parameter ";

	final StringBuilder sb = new StringBuilder(sequence);
	int nchunks = sequence.length() / width;
	// add up inserted new line chars
	nchunks = (nchunks + sequence.length()) / width;
	int nlineCharcounter = 0;
	for (int i = 1; i <= nchunks; i++) {
	    int insPos = width * i + nlineCharcounter;
	    // to prevent inserting new line in the very end of a sequence then
	    // it would have failed.
	    // Also covers the case when the sequences shorter than width
	    if (sb.length() <= insPos) {
		break;
	    }
	    sb.insert(insPos, "\n");
	    nlineCharcounter++;
	}
	return sb.toString();
    }

    /**
     * 
     * @return sequence length
     */
    public int getLength() {
	return this.sequence.length();
    }

    /**
     * Same as oneLineFasta
     */
    @Override
    public String toString() {
	return this.getOnelineFasta();
    }

    @Override
    public boolean equals(Object obj) {
	if (obj == null) {
	    return false;
	}
	if (!(obj instanceof FastaSequence)) {
	    return false;
	}
	FastaSequence fs = (FastaSequence) obj;
	if (!fs.getId().equals(this.getId())) {
	    return false;
	}
	if (!fs.getSequence().equalsIgnoreCase(this.getSequence())) {
	    return false;
	}
	return true;
    }

}
