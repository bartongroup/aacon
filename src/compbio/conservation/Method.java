package compbio.conservation;

	enum Method { KABAT_SCORE, JORES_SCORE, SCHNEIDER_SCORE, SHENKIN_SCORE, GERSTEIN_SCORE, TAYLOR_SCORE_GAPS, TAYLOR_SCORE_NO_GAPS, ZVELIBIL_SCORE, KARLIN_SCORE, ARMON_SCORE, THOMPSON_SCORE, NOT_LANCET_SCORE, MIRNY_SCORE, WILLIAMSON_SCORE, LANDGRAF_SCORE, SANDER_SCORE, VALDAR_SCORE;

	static Method getMethod(String meth) {

		meth = meth.trim().toLowerCase();

		if(meth.equalsIgnoreCase(KABAT_SCORE.toString())) {
	
			return KABAT_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(JORES_SCORE.toString())) {
		
			return JORES_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(SCHNEIDER_SCORE.toString())) {
		
			return SCHNEIDER_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(SHENKIN_SCORE.toString())) {
		
			return SHENKIN_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(GERSTEIN_SCORE.toString())) {
		
			return GERSTEIN_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(TAYLOR_SCORE_GAPS.toString())) {
		
			return TAYLOR_SCORE_GAPS;
		
		}
	
		if(meth.equalsIgnoreCase(TAYLOR_SCORE_NO_GAPS.toString())) {
		
			return TAYLOR_SCORE_NO_GAPS;
		
		}
	
		if(meth.equalsIgnoreCase(ZVELIBIL_SCORE.toString())) {
		
		return ZVELIBIL_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(KARLIN_SCORE.toString())) {
		
			return KARLIN_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(ARMON_SCORE.toString())) {
		
			return ARMON_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(THOMPSON_SCORE.toString())) {
		
			return THOMPSON_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(NOT_LANCET_SCORE.toString())) {
		
			return NOT_LANCET_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(MIRNY_SCORE.toString())) {
		
			return MIRNY_SCORE;
		
		}
		
		if(meth.equalsIgnoreCase(WILLIAMSON_SCORE.toString())) {
			
			return WILLIAMSON_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(LANDGRAF_SCORE.toString())) {
		
			return LANDGRAF_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(SANDER_SCORE.toString())) {
		
			return SANDER_SCORE;
		
		}
	
		if(meth.equalsIgnoreCase(VALDAR_SCORE.toString())) {
		
			return VALDAR_SCORE;
		
		}
		
	return null;

	}
	
	static void supportedMethods() {
		
		System.out.println("Supported methods:");
		
		for(Method method : Method.values()) {
			
			System.out.println(method.toString());
		}
	}

}
