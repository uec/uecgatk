package org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel;

import org.broadinstitute.sting.utils.exceptions.UserException;

/*
 * Bis-SNP/BisSNP: It is a genotyping and methylation calling in bisulfite treated 
 * massively parallel sequencing (Bisulfite-seq and NOMe-seq) on Illumina platform
 * Copyright (C) <2011>  <Yaping Liu: lyping1986@gmail.com>

 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

public class BisulfiteUserException extends UserException {

	public BisulfiteUserException(String msg) {
		super(msg);
		// TODO Auto-generated constructor stub
	}

	public BisulfiteUserException(String msg, Throwable e) {
		super(msg, e);
		// TODO Auto-generated constructor stub
	}
	
	

}
