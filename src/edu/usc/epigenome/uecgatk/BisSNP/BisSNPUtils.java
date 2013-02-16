package edu.usc.epigenome.uecgatk.BisSNP;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
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
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
	
	public static boolean isCytosine(byte base, boolean bisulfiteConversionSpace, boolean secondPair)
	{
		char refC = (char) base;
		boolean out;
		
		if (bisulfiteConversionSpace)
		{
			if(secondPair){
				out = ((refC == 'G') || (refC == 'A'));
			}
			else{
				out = ((refC == 'C') || (refC == 'T'));
			}
			
		}
		else
		{
			if(secondPair){
				out = (refC == 'G');
			}
			else{
				out = (refC == 'C');
			}
			
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
	
	public static boolean isCytosine(int pos, String seqStr, boolean bisulfiteConversionSpace, boolean secondPair)
	{
		char refC = seqStr.charAt(pos);
		
		boolean out;
		
		if (bisulfiteConversionSpace)
		{
			if(secondPair){
				out = ((refC == 'G') || (refC == 'A'));
			}
			else{
				out = ((refC == 'C') || (refC == 'T'));
			}
			
		}
		else
		{
			if(secondPair){
				out = (refC == 'G');
			}
			else{
				out = (refC == 'C');
			}
			
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
	
	public static boolean isHomoC(VariantContext vc) {
  	   	 if(vc.hasGenotypes()){
  	   		Iterator<Genotype> it = vc.getGenotypes().iterator();
  	   		 while(it.hasNext()){
  	   			Genotype tmp = it.next();
  	   			//System.err.println(tmp.getGenotypeString());
  	   			if(tmp.getGenotypeString().equalsIgnoreCase("C/C") || tmp.getGenotypeString().equalsIgnoreCase("G/G")){
  	   				return true;
  	   			}
  	   		 }
  	     }
  	     return false;
  	}
	
	public static boolean isHetCpattern(VariantContext vc) {  //het at C sites or G sites.
 	   	 if(vc.hasGenotypes()){
 	   		Iterator<Genotype> it = vc.getGenotypes().iterator();
 	   		 while(it.hasNext()){
 	   			Genotype tmp = it.next();
 	   			//System.err.println(tmp.getGenotypeString());
 	   			if(tmp.getGenotypeString().equalsIgnoreCase("C/C") || tmp.getGenotypeString().equalsIgnoreCase("G/G")){
 	   				return true;
 	   			}
 	   		 }
 	     }
 	     return false;
 	}
	/*
	public static boolean isHetCpg_at_C(VariantContext vc) {
	   	 if(vc.hasGenotypes()){
	   		Iterator<Genotype> it = vc.getGenotypes().iterator();
	   		 while(it.hasNext()){
	   			Genotype tmp = it.next();
	   			if(tmp.isHom())
	   				return false;
	   			//System.err.println(tmp.getGenotypeString());
	   			if(tmp.hasAttribute(BisulfiteVCFConstants.BEST_C_PATTERN)){
	   				byte[] bases = tmp.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").getBytes();
	   				if(bases.length < 2)
	   					return false;
	   				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[0],BaseUtils.C) && BaseUtils.basesAreEqual(bases[1], BaseUtils.G)){
	   					return true;
	   				}
	   				
	   			}
	   		 }
	     }
	     return false;
	}
	*/
	/*
	public static boolean isHetCpg_at_G(VariantContext vc) {
	   	 if(vc.hasGenotypes()){
	   		Iterator<Genotype> it = vc.getGenotypes().iterator();
	   		 while(it.hasNext()){
	   			Genotype tmp = it.next();
	   			if(tmp.isHom()){
	   				if(vc.hasAttribute(BisulfiteVCFConstants.CYTOSINE_TYPE)){
	   					System.err.println(vc.getAttributeAsString(BisulfiteVCFConstants.CYTOSINE_TYPE, "."));
	   					if(vc.getAttributeAsString(BisulfiteVCFConstants.CYTOSINE_TYPE, ".").equalsIgnoreCase("[CG, CH, C]"))
	   						return true;	
	   				}
	   			}
	   		 }
	     }
	     return false;
	}
	*/
	
	public static boolean isHetCpg_at_C(VariantContext vc) {
 	   	 if(vc.hasGenotypes()){
 	   		Iterator<Genotype> it = vc.getGenotypes().iterator();
 	   		 while(it.hasNext()){
 	   			Genotype tmp = it.next();
 	   			if(tmp.isHom())
 	   				return false;
 	   			//System.err.println(tmp.getGenotypeString());
 	   			if(tmp.hasAttribute(BisulfiteVCFConstants.BEST_C_PATTERN)){
 	   				byte[] bases = tmp.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").getBytes();
 	   				if(bases.length < 2)
 	   					return false;
 	   				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[0],BaseUtils.C) && BaseUtils.basesAreEqual(bases[1], BaseUtils.G)){
 	   					return true;
 	   				}
 	   				
 	   			}
 	   		 }
 	     }
 	     return false;
 	}
	
	public static boolean isHetCpg_at_G(VariantContext vc) {
	   	 if(vc.hasGenotypes()){
	   		Iterator<Genotype> it = vc.getGenotypes().iterator();
	   		 while(it.hasNext()){
	   			Genotype tmp = it.next();
	   			if(tmp.isHom()){
	   				if(tmp.hasAttribute(BisulfiteVCFConstants.BEST_C_PATTERN)){
	 	   				byte[] bases = tmp.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").getBytes();
	 	   			if(bases.length < 2)
 	   					return false;
	 	   				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[1],BaseUtils.G) && BaseUtils.basesAreEqual(bases[0], BaseUtils.C)){
	 	   					return true;
	 	   				}
	 	   				
	 	   			}
	   			}
	   		 }
	     }
	     return false;
	}
	
	public static boolean isHetCph_at_C(VariantContext vc) {
	   	 if(vc.hasGenotypes()){
	   		Iterator<Genotype> it = vc.getGenotypes().iterator();
	   		 while(it.hasNext()){
	   			Genotype tmp = it.next();
	   			if(tmp.isHom())
	   				return false;
	   			//System.err.println(tmp.getGenotypeString());
	   			if(tmp.hasAttribute(BisulfiteVCFConstants.BEST_C_PATTERN)){
	   				byte[] bases = tmp.getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").getBytes();
	   				if(bases.length < 2)
 	   					return false;
	   				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[0],BaseUtils.C) && bases[1]==BaseUtilsMore.H){
	   					return true;
	   				}
	   				
	   			}
	   		 }
	     }
	     return false;
	}
	
	
	public static boolean isRefCytosinePattern(ReferenceContext refContext, String cytosinePattern){
		if(isRefCytosinePattern(refContext, cytosinePattern, false) || isRefCytosinePattern(refContext, cytosinePattern, true))
			return true;
		return false;
	}
	
	public static boolean isRefCytosinePattern(ReferenceContext refContext, String cytosinePattern, boolean negStrand){
		
		byte[] bases = cytosinePattern.getBytes();
		int matches=0;
		if(negStrand){
			byte[] refBytesRev = new byte[cytosinePattern.length()];
			
			for(int i=0; i < cytosinePattern.length(); i++){
				GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() - i );
				if( !refContext.getWindow().containsP(loc) ){
					refBytesRev[i]=BaseUtils.NO_CALL_INDEX;
					continue;
				}

				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),loc, refContext.getWindow(),refContext.getBases());
				refBytesRev[i] = tmpRef.getBase();
			}
			refBytesRev = BaseUtilsMore.simpleReverse(refBytesRev);
			for(int i=0; i<bases.length; i++){
				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[i], refBytesRev[i])){
					matches++;
				}
			}
			if(matches == bases.length)
				return true;
		}
		else{
			byte[] refBytesFwd = new byte[cytosinePattern.length()];
			for(int i=0; i < cytosinePattern.length(); i++){
				GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() + i );
				if( !refContext.getWindow().containsP(loc) ){
					refBytesFwd[i]=BaseUtils.NO_CALL_INDEX;
					continue;
				}
					
				
				ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),loc, refContext.getWindow(),refContext.getBases());
				refBytesFwd[i] = tmpRef.getBase();
			}
			
			for(int i=0; i<bases.length; i++){
				if(BaseUtilsMore.iupacCodeEqualNotConsiderMethyStatus(bases[i], refBytesFwd[i])){
					matches++;
				}
			}
			if(matches == bases.length)
				return true;
		}
		
		return false;
	}
	
	public static boolean isRefCpg(ReferenceContext refContext){
		byte[] refBytes = new byte[3];
		for(int i=-1, index=0; i <= 1; i++, index++){
			GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() + i );
			if( !refContext.getWindow().containsP(loc) ){
				refBytes[index]=BaseUtils.NO_CALL_INDEX;
				continue;
			}
				
			
			ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),loc, refContext.getWindow(),refContext.getBases());
			refBytes[index] = tmpRef.getBase();
		}
		if(BaseUtils.basesAreEqual(refBytes[0], BaseUtils.C) && BaseUtils.basesAreEqual(refBytes[1], BaseUtils.G)){
			return true;
		}	
		else if(BaseUtils.basesAreEqual(refBytes[1], BaseUtils.C) && BaseUtils.basesAreEqual(refBytes[2], BaseUtils.G)){
			return true;
		}
		else{
			return false;
		}
	}
	
	public static boolean isRefCph(ReferenceContext refContext){
		byte[] refBytes = new byte[3];
		for(int i=-1, index=0; i <= 1; i++, index++){
			GenomeLoc loc = refContext.getGenomeLocParser().createGenomeLoc(refContext.getLocus().getContig(), refContext.getLocus().getStart() + i );
			if( !refContext.getWindow().containsP(loc) ){
				refBytes[index]=BaseUtils.NO_CALL_INDEX;
				continue;
			}
				
			
			ReferenceContext tmpRef = new ReferenceContext(refContext.getGenomeLocParser(),loc, refContext.getWindow(),refContext.getBases());
			refBytes[index] = tmpRef.getBase();
		}
		if(BaseUtils.basesAreEqual(refBytes[1], BaseUtils.G) && (!BaseUtils.basesAreEqual(refBytes[0], BaseUtils.C) && !BaseUtils.basesAreEqual(refBytes[0], BaseUtils.N) && !BaseUtils.basesAreEqual(refBytes[0], BaseUtils.NO_CALL_INDEX))){
			//System.err.println(refContext.getLocus().getStart() + "\tneg:" + new String(refBytes));
			return true;
		}	
		else if(BaseUtils.basesAreEqual(refBytes[1], BaseUtils.C) && (!BaseUtils.basesAreEqual(refBytes[2], BaseUtils.G) && !BaseUtils.basesAreEqual(refBytes[2], BaseUtils.N) && !BaseUtils.basesAreEqual(refBytes[2], BaseUtils.NO_CALL_INDEX))){
			//System.err.println(refContext.getLocus().getStart() + "\tpos:" + new String(refBytes));
			return true;
		}
		else{
			//System.err.println(refContext.getLocus().getStart() + "\tnot:" + new String(refBytes));
			return false;
		}
	}
}
