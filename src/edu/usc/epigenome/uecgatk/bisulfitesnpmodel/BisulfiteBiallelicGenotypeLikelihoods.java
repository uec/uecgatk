package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import org.broad.tribble.util.variantcontext.Allele;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;

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

public class BisulfiteBiallelicGenotypeLikelihoods extends
		BiallelicGenotypeLikelihoods {

	private Integer[] CytosineStatus;
	
	public BisulfiteBiallelicGenotypeLikelihoods(String sample, Allele A,
			Allele B, double log10aaLikelihoods, double log10abLikelihoods,
			double log10bbLikelihoods, int depth) {
		super(sample, A, B, log10aaLikelihoods, log10abLikelihoods,
				log10bbLikelihoods, depth);
		// TODO Auto-generated constructor stub
	}
	
	public BisulfiteBiallelicGenotypeLikelihoods(String sample, Allele A,
			Allele B, double log10aaLikelihoods, double log10abLikelihoods,
			double log10bbLikelihoods, int depth, Integer[] CytosineStatus) {
		super(sample, A, B, log10aaLikelihoods, log10abLikelihoods,
				log10bbLikelihoods, depth);
		this.CytosineStatus = CytosineStatus;
		// TODO Auto-generated constructor stub
	}
	
	public Integer[] getCytosineStatus() {
        return this.CytosineStatus;
    }

}
