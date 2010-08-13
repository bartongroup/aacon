package compbio.conservation;

/**
 * Thrown when there is an illegal(not an amino acid or known gap character) in the
 * sequence list or character array fed into AminoAcidMatrix constructor)
 * @author agolicz
 *
 */

	class NotAnAminoAcidException extends RuntimeException {
	
	NotAnAminoAcidException(String message) {
		
		super(message);
	}

}
