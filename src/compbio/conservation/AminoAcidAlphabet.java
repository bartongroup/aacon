package compbio.conservation;

public class AminoAcidAlphabet {

	private char[] alph;

	public AminoAcidAlphabet() {

	alph = new char[21];

	alph[0] = 'R';
	alph[1] = 'H';
	alph[2] = 'K';
	alph[3] = 'D';
	alph[4] = 'E';
	alph[5] = 'S';
	alph[6] = 'T';
	alph[7] = 'N';
	alph[8] = 'Q';
	alph[9] = 'C';
	alph[10] = 'G';
	alph[11] = 'P';
	alph[12] = 'A';
	alph[13] = 'I';
	alph[14] = 'L';
	alph[15] = 'M';
	alph[16] = 'F';
	alph[17] = 'W';
	alph[18] = 'Y';
	alph[19] = 'V';
	alph[20] = '-';

	}

	public char[] getAlphabet() {

	return alph;

	}
	
	public int length() {
		
	int len = 21;
	
	return len;
	
	}

}
