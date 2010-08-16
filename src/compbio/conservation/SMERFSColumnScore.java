package compbio.conservation;

/**
 * Enumeration listing two possibilities of how to give window scores to the
 * columns while calculating SMERFS.
 * 
 * MAX_SCORE gives the highest core of all the windows the column belongs to
 * MID_SCORE give the window score to the column in the middle
 * 
 * @author agolicz
 * 
 */

public enum SMERFSColumnScore {
    MAX_SCORE, MID_SCORE;

    static SMERFSColumnScore getSMERFSColumnScore(String score) {

	score = score.trim().toLowerCase();

	if (score.equalsIgnoreCase(SMERFSColumnScore.MAX_SCORE.toString())) {

	    return SMERFSColumnScore.MAX_SCORE;
	}

	if (score.equalsIgnoreCase(SMERFSColumnScore.MID_SCORE.toString())) {

	    return SMERFSColumnScore.MID_SCORE;
	}

	return null;
    }

    static void supportedSMERFSColumnSores() {

	System.out.println("Supported formats:");

	for (SMERFSColumnScore score : SMERFSColumnScore.values()) {

	    System.out.println(score.toString());
	}
    }

}
