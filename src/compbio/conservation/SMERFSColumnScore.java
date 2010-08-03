package compbio.conservation;

public enum SMERFSColumnScore { MAX_SCORE, MID_SCORE;

	static SMERFSColumnScore getSMERFSColumnScore(String score) {
	
		score = score.trim().toLowerCase();
		
		if(score.equalsIgnoreCase(SMERFSColumnScore.MAX_SCORE.toString())) {
		
			return SMERFSColumnScore.MAX_SCORE;
		}
	
		if(score.equalsIgnoreCase(SMERFSColumnScore.MID_SCORE.toString())) {
		
			return SMERFSColumnScore.MID_SCORE;
		}
	
		return null;
	}

	static void supportedSMERFSColumnSores() {
	
		System.out.println("Supported formats:");
	
		for(SMERFSColumnScore score : SMERFSColumnScore.values()) {
		
			System.out.println(score.toString());
		}
	}


}
