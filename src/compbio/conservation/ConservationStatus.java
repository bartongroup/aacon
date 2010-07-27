package compbio.conservation;

public enum ConservationStatus { IDENTICAL, CONSERVED, NOT_CONSERVED;
	
	static ConservationStatus getStatus(double score, Method method) {
	
	if (score == 1) {
		
		return ConservationStatus.IDENTICAL;
	}
	
	else {
		
		if (score >= SubFamiliesConservation.getTreshhold(method)) {
			
			return ConservationStatus.CONSERVED;
			
		}

		else {
			
			return ConservationStatus.NOT_CONSERVED;
		}

	}
	
	}
	
	static String stringReps(ConservationStatus con) {
		
		if (con == ConservationStatus.IDENTICAL) {
			
			return IDENTICAL.toString();
		}
		
		if (con == CONSERVED) {
			
			return CONSERVED.toString();
		}
		
		if (con == NOT_CONSERVED) {
			
			return NOT_CONSERVED.toString();
		}
		
		return null;
	}


}
