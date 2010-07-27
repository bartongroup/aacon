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

}
