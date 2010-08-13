package compbio.conservation;

/**
 * Thrown when an argument provided as a gap character can not be parsed as char.
 * @author agolicz
 *
 */

public class IllegalGapCharacterException extends IllegalArgumentException{
	
	public IllegalGapCharacterException(String message) {
		
		super(message);
		
	}

}
