package compbio.conservation;

import java.io.FileNotFoundException;
import java.util.concurrent.Callable;

import compbio.util.Timer;

/**
 * TODO to complete
 * 
 * @author pvtroshin
 */
public class MethodWrapper implements Callable<MethodWrapper> {

	double[] conservation = null;
	final Method method;
	private final ConservationScores2 scores;
	private final boolean normalize;
	private int SMERFSWidth = SMERFSColumnScore.DEFAULT_WINDOW_SIZE;
	private double SMERFSGapTreshold = SMERFSColumnScore.DEFAULT_GAP_THRESHOLD;
	private SMERFSColumnScore colScoreSchema = SMERFSColumnScore.MID_SCORE;

	final Timer timer;

	public MethodWrapper(Method method, ConservationScores2 scores,
			boolean normalize, Timer timer) throws FileNotFoundException {
		this.method = method;
		this.scores = scores;
		this.normalize = normalize;
		this.timer = timer;
	}

	private boolean isCustomSmerfsParams() {
		if (this.colScoreSchema == SMERFSColumnScore.MID_SCORE
				&& this.SMERFSGapTreshold == SMERFSColumnScore.DEFAULT_GAP_THRESHOLD
				&& this.SMERFSWidth == SMERFSColumnScore.DEFAULT_WINDOW_SIZE) {
			return false;
		}
		return true;
	}

	public MethodWrapper(ConservationScores2 scores, boolean normalize,
			int SMERFSWidth, SMERFSColumnScore colScoreSchema,
			double SMERFSGapTreshold, Timer timer) throws FileNotFoundException {
		this(Method.SMERFS, scores, normalize, timer);
		this.colScoreSchema = colScoreSchema;
		this.SMERFSWidth = SMERFSWidth;
		this.SMERFSGapTreshold = SMERFSGapTreshold;
	}

	@Override
	public MethodWrapper call() throws Exception {
		timer.getStepTime();
		if (method == Method.SMERFS && isCustomSmerfsParams()) {
			this.conservation = ConservationClient.getSMERFS(scores
					.getAlignment(), this.SMERFSWidth, this.colScoreSchema,
					this.SMERFSGapTreshold, this.normalize);
		} else {
			this.conservation = scores.calculateScore(method, normalize);
		}
		timer.println(method.toString() + " " + timer.getStepTime() + " ms");
		return this;
	}
}
