package compbio.conservation;

	enum Method { KABAT, JORES, SCHNEIDER, SHENKIN, GERSTEIN, TAYLOR_GAPS, TAYLOR_NO_GAPS, ZVELIBIL, KARLIN, ARMON, THOMPSON, NOT_LANCET, MIRNY, WILLIAMSON, LANDGRAF, SANDER, VALDAR, SMERFS;

	static Method getMethod(String meth) {

		meth = meth.trim().toLowerCase();

		if(meth.equalsIgnoreCase(KABAT.toString())) {
	
			return KABAT;
		
		}
	
		if(meth.equalsIgnoreCase(JORES.toString())) {
		
			return JORES;
		
		}
	
		if(meth.equalsIgnoreCase(SCHNEIDER.toString())) {
		
			return SCHNEIDER;
		
		}
	
		if(meth.equalsIgnoreCase(SHENKIN.toString())) {
		
			return SHENKIN;
		
		}
	
		if(meth.equalsIgnoreCase(GERSTEIN.toString())) {
		
			return GERSTEIN;
		
		}
	
		if(meth.equalsIgnoreCase(TAYLOR_GAPS.toString())) {
		
			return TAYLOR_GAPS;
		
		}
	
		if(meth.equalsIgnoreCase(TAYLOR_NO_GAPS.toString())) {
		
			return TAYLOR_NO_GAPS;
		
		}
	
		if(meth.equalsIgnoreCase(ZVELIBIL.toString())) {
		
		return ZVELIBIL;
		
		}
	
		if(meth.equalsIgnoreCase(KARLIN.toString())) {
		
			return KARLIN;
		
		}
	
		if(meth.equalsIgnoreCase(ARMON.toString())) {
		
			return ARMON;
		
		}
	
		if(meth.equalsIgnoreCase(THOMPSON.toString())) {
		
			return THOMPSON;
		
		}
	
		if(meth.equalsIgnoreCase(NOT_LANCET.toString())) {
		
			return NOT_LANCET;
		
		}
	
		if(meth.equalsIgnoreCase(MIRNY.toString())) {
		
			return MIRNY;
		
		}
		
		if(meth.equalsIgnoreCase(WILLIAMSON.toString())) {
			
			return WILLIAMSON;
		
		}
	
		if(meth.equalsIgnoreCase(LANDGRAF.toString())) {
		
			return LANDGRAF;
		
		}
	
		if(meth.equalsIgnoreCase(SANDER.toString())) {
		
			return SANDER;
		
		}
	
		if(meth.equalsIgnoreCase(VALDAR.toString())) {
		
			return VALDAR;
		
		}
		
		if(meth.equalsIgnoreCase(SMERFS.toString())) {
			
			return SMERFS;
			
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
