package compbio.conservation;

import java.util.HashMap;
import java.util.Map;

class AminoAcidAlphabet {

    private final char[] alph;
    private final Map<Character, Integer> charCount;

    AminoAcidAlphabet() {
	charCount = new HashMap<Character, Integer>();

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

    Map<Character, Integer> calculateOccurance(char[] column) {
	for (char ch : column) {
	    addToOccurance(ch);
	}
	return charCount;
    }

    void addToOccurance(char ch) {
	Integer count = charCount.get(ch);
	if (count == null) {
	    charCount.put(ch, new Integer(1));
	} else {
	    charCount.put(ch, count + 1);
	}
    }

    public char[] getAlphabet() {

	return alph;

    }

    public int length() {

	int len = 21;

	return len;

    }

}
