package compbio.conservation;

import java.io.FileNotFoundException;
import java.util.concurrent.Callable;

import compbio.util.Timer;

/**
 * Wrapper for AA Conservation calculation methods and their results to enable
 * parallel method execution
 * 
 * @author pvtroshin
 */
public class MethodWrapper implements Callable<MethodWrapper> {

	double[] conservation = null;
	final Method method;
	private final Conservation scores;
	private final boolean normalize;
	private int SMERFSWidth = SMERFSColumnScore.DEFAULT_WINDOW_SIZE;
	private double SMERFSGapTreshold = SMERFSColumnScore.DEFAULT_GAP_THRESHOLD;
	private SMERFSColumnScore colScoreSchema = SMERFSColumnScore.MID_SCORE;

	final Timer timer;

	public MethodWrapper(Method method, Conservation scores,
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

	public MethodWrapper(Conservation scores, boolean normalize,
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
			Conservation conserv = new Conservation(scores
					.getAlignment(), normalize);
			this.conservation = conserv.getSMERFS(this.SMERFSWidth,
					this.colScoreSchema, this.SMERFSGapTreshold);
		} else {
			this.conservation = scores.calculateScore(method);
		}
		timer.println(method.toString() + " " + timer.getStepTime() + " ms");
		return this;
	}
}
