/*
 * Copyright (c) 2010 Agnieszka Golicz & Peter Troshin 
 * 
 * Amino Acid Conservation @version: 1.0 
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the Apache License version 2 as published by the
 * Apache Software Foundation This library is distributed in the hope that it
 * will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache
 * License for more details. A copy of the license is in apache_license.txt. It
 * is also available here: http://www.apache.org/licenses/LICENSE-2.0.txt 
 * Any republication or derived work distributed in source code form must 
 * include this copyright and license notice.
 * 
 */
package compbio.conservation;

import java.util.concurrent.Callable;

import compbio.data.sequence.Method;
import compbio.util.Timer;

/**
 * Wrapper for AA Conservation calculation methods and their results to enable
 * parallel method execution
 * 
 * @author Peter Troshin
 */
final class MethodWrapper implements Callable<MethodWrapper> {

	double[] conservation = null;
	final Method method;
	private final Conservation scores;

	final Timer timer;

	MethodWrapper(Method method, Conservation scores, Timer timer) {
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
