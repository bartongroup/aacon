package compbio.conservation;

import compbio.util.FastaSequence;

public class _SequenceToMatrix {

	FastaSequence fsequence; 
	
	public _SequenceToMatrix(FastaSequence fsequence) {
		this.fsequence = fsequence;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((fsequence == null) ? 0 : fsequence.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		_SequenceToMatrix other = (_SequenceToMatrix) obj;
		if (fsequence == null) {
			if (other.fsequence != null)
				return false;
		} else if (!fsequence.equals(other.fsequence))
			return false;
		return true;
	}
	
}
