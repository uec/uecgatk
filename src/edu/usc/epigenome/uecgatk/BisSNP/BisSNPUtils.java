package edu.usc.epigenome.uecgatk.BisSNP;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;

import edu.usc.epigenome.uecgatk.YapingWalker.NDRargumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteGenotyperEngine.OUTPUT_MODE;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BisulfiteSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel;

import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.BadBaseFilterBisulfite;

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

/*
 * provide easy access for some function inside BisSNP.
 *  e.g. judge if it is a type of cytosine pattern(HCG-2, or GCH-1) from the given pileup and corresponding reference seq, dbSNP. 
 * should provide methylation value, otherwise will use flat methylation value; and also likelihood ratio criteria, bisulfite conversion rate
 */

public class BisSNPUtils {
	
	private double FLAT_METHY_STATUS = 0.5;
	private NDRargumentCollection NAC;
	private BisulfiteDiploidSNPGenotypePriors genotypePriors;
	private BisulfiteArgumentCollection BAC;
	
	public BisSNPUtils(BisulfiteArgumentCollection BAC){
		this.NAC = new NDRargumentCollection();
		this.genotypePriors = new BisulfiteDiploidSNPGenotypePriors();
	}
	
	public BisSNPUtils(NDRargumentCollection NAC){
		this.NAC = NAC;
		this.genotypePriors = new BisulfiteDiploidSNPGenotypePriors();
	}
	
	
	
	public static String getRefGenomeVersion(){
		return BisulfiteGenotyper.getBAC().referenceGenomeVer;
	}
	
	public class methyStatus{
		DiploidGenotype genotype;
		double ratio;
		methyStatus(){
			
		}
	}
	
	public static boolean isCytosine(byte base, boolean bisulfiteConversionSpace)
	{
		char refC = (char) base;
		boolean out;
		
		if (bisulfiteConversionSpace)
		{
			out = ((refC == 'C') || (refC == 'T'));
		}
		else
		{
			out = (refC == 'C');
		}
		
		return out; 
	}
	
	public static boolean isCytosine(int pos, String seqStr, boolean bisulfiteConversionSpace)
	{
		char refC = seqStr.charAt(pos);
		
		boolean out;
		
		if (bisulfiteConversionSpace)
		{
			out = ((refC == 'C') || (refC == 'T'));
		}
		else
		{
			out = (refC == 'C');
		}
		
		return out; 
	}
	
	public static boolean isThymine(int pos, String seqStr)
	{
		char seqC = seqStr.charAt(pos);
		
		return (seqC == 'T'); 
	}
	
	public static boolean isThymineInCytosinePos(int pos, String seqStr, byte refBase)
	{
		char seqC = seqStr.charAt(pos);
		
		return ((seqC == 'T') && (refBase == 'C')); 
	}
	
	
}
