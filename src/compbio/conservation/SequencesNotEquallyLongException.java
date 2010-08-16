package compbio.conservation;

/**
 * Thrown when not all the sequences in the fasta file are of the same length.
 * 
 * @author agolicz
 * 
 */

public class SequencesNotEquallyLongException extends RuntimeException {

    SequencesNotEquallyLongException(String message) {

	super(message);
    }

}
