package compbio.conservation;

public enum Format { RESULT_WITH_ALIGNMENT, RESULT_NO_ALIGNMENT;

	static Format getFormat(String format) {
		
		format = format.trim().toLowerCase();
		
		if(format.equalsIgnoreCase(Format.RESULT_WITH_ALIGNMENT.toString())) {
			
			return Format.RESULT_WITH_ALIGNMENT;
		}
		
		if(format.equalsIgnoreCase(Format.RESULT_NO_ALIGNMENT.toString())) {
			
			return Format.RESULT_NO_ALIGNMENT;
		}
		
		return null;
	}
	
	static void supportedFormats() {
		
		System.out.println("Supported formats:");
		
		for(Format format : Format.values()) {
			
			System.out.println(format.toString());
		}
	}
	
}
