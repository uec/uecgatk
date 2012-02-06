package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.util.Map;

import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext.StratifiedContextType;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.genotyper.BiallelicGenotypeLikelihoods;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;

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

public class NonRefDependSNPGenotypeLikelihoodsCalculationModel extends
		GenotypeLikelihoodsCalculationModel {

	public enum MethylSNPModel {
        BM,
        GM,
        NM
    }
	
	public NonRefDependSNPGenotypeLikelihoodsCalculationModel(
			UnifiedArgumentCollection UAC, Logger logger) {
		super(UAC, logger);
		// TODO Auto-generated constructor stub
	}

	@Override
	public Allele getLikelihoods(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, StratifiedAlignmentContext> contexts,
			StratifiedContextType contextType, GenotypePriors priors,
			Map<String, BiallelicGenotypeLikelihoods> GLs,
			Allele alternateAlleleToUse) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public Allele getBsLikelihoods(RefMetaDataTracker tracker,
			ReferenceContext ref,
			Map<String, StratifiedAlignmentContext> contexts,
			StratifiedContextType contextType,
			Map<String, BisulfiteBiallelicGenotypeLikelihoods> GLs,
			Allele alternateAlleleToUse) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void initialize(CytosineTypeStatus cts, BisulfiteArgumentCollection BAC, boolean autoEstimateC, boolean secondIteration){

	}

}
