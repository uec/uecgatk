/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 20, 2012 11:46:27 AM
 * 
 */
public class BisulfiteVariantCallContext {

	private HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs;
	private VariantContext vc;
	private summaryAcrossReadsGroup summary;
	
	public boolean isC = false;
	public boolean isCpg = false;
	public boolean isHetCpg = false;
	public AlignmentContext rawContext = null;
	public ReferenceContext ref = null;
    public boolean confidentlyCalled = false;
    // Should this site be emitted?
    public boolean shouldEmit = true;
	/**
	 * 
	 */
	public BisulfiteVariantCallContext(HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs, VariantContext vc, AlignmentContext rawContext, ReferenceContext ref) {
		// TODO Auto-generated constructor stub
		this.rawContext = rawContext;
		this.ref = ref;
		this.BCGLs = BCGLs;
		this.vc = vc;
		setSummaryAcrossReadsGroup();
		
	}
	
	public HashMap<String, BisulfiteContextsGenotypeLikelihoods> getBisulfiteContextsGenotypeLikelihoods(){
		return BCGLs;
	}
	
	public VariantContext getVariantContext(){
		return vc;
	}
	
	public summaryAcrossReadsGroup getSummaryAcrossRG(){
		return summary;
	}
	
	public boolean isHetSnp() {
   	 if(this.vc.hasGenotypes()){
   		Iterator<Genotype> it = vc.getGenotypes().iterator();
   		 while(it.hasNext()){
   			 if(it.next().isHet())
   				 return true;
   		 }
     }
        return false;
   }
	
	public boolean isVariant() {
	   	 if(this.vc.hasGenotypes()){
	   		Iterator<Genotype> it = vc.getGenotypes().iterator();
	   		 while(it.hasNext()){
	   			Genotype tmp = it.next();
	   			 if(tmp.isHet() || tmp.isHomVar())
	   				 return true;
	   		 }
	     }
	        return false;
	   }
	
	private void setSummaryAcrossReadsGroup(){
		summary = new summaryAcrossReadsGroup();
//		boolean first = true;
		summary.cytosinePatternConfirmedSet = new HashSet<String>();
		if(!BCGLs.isEmpty()){
			for(BisulfiteContextsGenotypeLikelihoods BCGL : BCGLs.values()){
				summary.numC += BCGL.getNumOfCReadsInBisulfiteCStrand();
				summary.numT += BCGL.getNumOfTReadsInBisulfiteCStrand();
				summary.numA += BCGL.getNumOfAReadsInGenotypeGStrand();
				summary.numG += BCGL.getNumOfGReadsInGenotypeGStrand();
				for(String cytosinePattern : BCGL.getCytosineParameters().keySet()){
					
					if(BCGL.getCytosineParameters().get(cytosinePattern).isCytosinePattern){
						summary.cytosinePatternStrand = BCGL.getCytosineParameters().get(cytosinePattern).cytosineStrand;
						
			            	
						summary.cytosinePatternConfirmedSet.add(cytosinePattern);
						isC = true;
						if(cytosinePattern.equalsIgnoreCase("CG")){
							isCpg = true;
							if(BCGL.getCytosineParameters().get("CG").isHeterozygousCytosinePattern)
								isHetCpg = true;
						}
							
						
					}
					
				}
				
			}
			summary.cytosinePatternConfirmedList = summary.cytosinePatternConfirmedSet.toString();
		}
		
	}
	
	public class summaryAcrossReadsGroup{
		public int numC = 0;
		public int numT = 0;
		public int numG = 0;
		public int numA = 0;
		public char cytosinePatternStrand;
		public String cytosinePatternConfirmedList;
		public HashSet<String> cytosinePatternConfirmedSet;
	}

}
