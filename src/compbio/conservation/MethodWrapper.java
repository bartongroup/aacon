package compbio.conservation;

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

	final Timer timer;

	public MethodWrapper(Method method, Conservation scores, Timer timer) {
		this.method = method;
		this.scores = scores;
		this.timer = new Timer(timer);
	}

	@Override
	public MethodWrapper call() throws Exception {
		assert method != Method.SMERFS : " Must use separate method to calculate "
				+ "SMERFS to avoid thread contantion";

		timer.getStepTime();
		this.conservation = scores.calculateScore(method);
		timer.println(method.toString() + " " + timer.getStepTime() + " ms");
		return this;
	}
}
