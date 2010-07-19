package compbio.conservation;

//class for objects that represent amino acid occurnce 
//in particulat columns

public class _AminoAcidOccurance {

    char id;

    int freq;

    public _AminoAcidOccurance() {

	id = ' ';

	freq = 0;

    }

    public _AminoAcidOccurance(char id) {
	this.id = id;
	freq = 0;
    }

    public void addToOccurance(int a) {
	assert a > 0;
	freq = freq + a;

    }

    public int getOccurance() {

	return freq;

    }

    public char getId() {

	return id;

    }

}
