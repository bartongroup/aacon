package compbio.conservation;

public enum Smerfs {
    SMERFS;

    static Smerfs getSmerfs(String smerfs) {

	smerfs = smerfs.trim().toLowerCase();

	if (smerfs.equalsIgnoreCase(Smerfs.SMERFS.toString())) {

	    return SMERFS;

	}

	return null;

    }

}
