package compbio.conservation;

	enum Method { kabatScore, joresScore, schneiderScore, shenkinScore, gersteinScore, SmallestTaylorSetGaps, SmallestTaylorSetNoGaps, zvelibilScore, karlinScore, armonScore, thompsonScore, notLancetScore, mirnyScore, williamsonScore, landgrafScore, sanderScore, valdarScore;

	static Method getMethod(String meth) {

		meth = meth.trim().toLowerCase();

		if(meth.equalsIgnoreCase(kabatScore.toString())) {
	
			return kabatScore;
		
		}
	
		if(meth.equalsIgnoreCase(joresScore.toString())) {
		
			return joresScore;
		
		}
	
		if(meth.equalsIgnoreCase(schneiderScore.toString())) {
		
			return schneiderScore;
		
		}
	
		if(meth.equalsIgnoreCase(shenkinScore.toString())) {
		
			return shenkinScore;
		
		}
	
		if(meth.equalsIgnoreCase(gersteinScore.toString())) {
		
			return gersteinScore;
		
		}
	
		if(meth.equalsIgnoreCase(SmallestTaylorSetGaps.toString())) {
		
			return SmallestTaylorSetGaps;
		
		}
	
		if(meth.equalsIgnoreCase(SmallestTaylorSetNoGaps.toString())) {
		
			return SmallestTaylorSetNoGaps;
		
		}
	
		if(meth.equalsIgnoreCase(zvelibilScore.toString())) {
		
		return zvelibilScore;
		
		}
	
		if(meth.equalsIgnoreCase(karlinScore.toString())) {
		
			return karlinScore;
		
		}
	
		if(meth.equalsIgnoreCase(armonScore.toString())) {
		
			return armonScore;
		
		}
	
		if(meth.equalsIgnoreCase(thompsonScore.toString())) {
		
			return thompsonScore;
		
		}
	
		if(meth.equalsIgnoreCase(notLancetScore.toString())) {
		
			return notLancetScore;
		
		}
	
		if(meth.equalsIgnoreCase(mirnyScore.toString())) {
		
			return mirnyScore;
		
		}
		
		if(meth.equalsIgnoreCase(williamsonScore.toString())) {
			
			return williamsonScore;
		
		}
	
		if(meth.equalsIgnoreCase(landgrafScore.toString())) {
		
			return landgrafScore;
		
		}
	
		if(meth.equalsIgnoreCase(sanderScore.toString())) {
		
			return sanderScore;
		
		}
	
		if(meth.equalsIgnoreCase(valdarScore.toString())) {
		
			return valdarScore;
		
		}
		
	return null;

	}

}
