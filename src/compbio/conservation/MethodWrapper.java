package compbio.conservation;

import java.util.concurrent.Callable;

/**
 * TODO to complete
 * 
 * @author pvtroshin
 * 
 */
public class MethodWrapper implements Callable<double[]> {

    private final String method;
    private final ConservationScores2 scores;
    private final boolean normalize;

    private final int SMERFSWidth;
    private final double SMERFSGapTreshold;

    public MethodWrapper(String method, ConservationScores2 scores,
	    boolean normalize) {
	this.method = method;
	this.scores = scores;
	this.normalize = normalize;
	this.SMERFSWidth = -1;
	this.SMERFSGapTreshold = -1;
    }

    public MethodWrapper(ConservationScores2 scores, boolean normalize,
	    int SMERFSWidth, double SMERFSGapTreshold) {
	this.method = Method.SMERFS.toString();
	this.scores = scores;
	this.normalize = normalize;
	this.SMERFSWidth = SMERFSWidth;
	this.SMERFSGapTreshold = SMERFSGapTreshold;
    }

    @Override
    public double[] call() throws Exception {
	double[] result = ParallelConservationClient.getMethod(method
		.toString(), scores, normalize);

	return result;
    }

}
