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
package compbio.common;

/**
 * Thrown when there is an illegal(not an amino acid or known gap character) in
 * the sequence list or character array fed into AminoAcidMatrix constructor)
 * 
 * @author Agnieszka Golicz
 */
public class NotAnAminoAcidException extends RuntimeException {

	public NotAnAminoAcidException(String message) {

		super(message);
	}
}
