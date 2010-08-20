package compbio.common;

/**
 * Thrown when there is an illegal(not an amino acid or known gap character) in
 * the sequence list or character array fed into AminoAcidMatrix constructor)
 * 
 * @author agolicz
 */
public class NotAnAminoAcidException extends RuntimeException {

	public NotAnAminoAcidException(String message) {

		super(message);
	}
}
