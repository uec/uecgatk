/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.picard.reference.IndexedFastaSequenceFile;

import org.broad.tribble.bed.BEDFeature;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.IntervalStratification;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.SortableJexlVCMatchExp;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.VariantEvalUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import edu.usc.epigenome.uecgatk.bisulfiteIndels.BisulfiteRealignerTargetCreator;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 5, 2012 5:09:11 PM
 * 
 */
public class vcfComparisonWalker extends RodWalker<Integer,Integer>{

	/**
	 * 
	 */
	
	/**
     * The variant file(s) to evaluate.
     */
    @Input(fullName="eval", shortName = "eval", doc="Input evaluation file(s)", required=true)
    public RodBinding<VariantContext> eval;

    /**
     * The variant file(s) to compare against.
     */
    @Input(fullName="comp", shortName = "comp", doc="Input comparison file(s)", required=false)
    public RodBinding<VariantContext> comp;
  //  private List<RodBinding<VariantContext>> comps = new ArrayList<RodBinding<VariantContext>>();

    @Argument(fullName="homCytosineComp", shortName="homC", doc="compare homozygous Cytosines, otherwise compare heterozygous SNP", required=false)
    protected boolean homC = false;

    
    @Argument(fullName="minQualityEval", shortName="mq", doc="Minimum genotyping quality for evaluation vcf file", required=false)
    protected double mq = 20.0;
    
    @Argument(fullName="minQualityStandard", shortName="mqComp", doc="Minimum genotyping quality for standard vcf file", required=false)
    protected double mqComp = 50.0;
    
 private int eval_num=0;
 private int eval_num_filtered=0;
 private int eval_num_no_call=0;
 private int comp_num=0;
 private int tp=0;
 private int fp=0;
 private int fn=0;
 private int tn=0;
    
    public void initialize() {
        
        // maintain the full list of comps
    //    comps.addAll(compsProvided);
        
        
    }

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Integer map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		// TODO Auto-generated method stub
		if(tracker==null)
			return null;
		List<VariantContext> comp_bindings = tracker.getValues(comp);
		if(!comp_bindings.isEmpty()){
			
			
			VariantContext vc_comp = comp_bindings.get(0); 
			comparison(vc_comp,tracker);
			
			comp_num++;
		}
		
		 
		
		
		
        return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public Integer reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}

	
	

	

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Integer reduce(Integer value, Integer sum) {
		// TODO Auto-generated method stub
		return null;
	}

	public void onTraversalDone(Integer result) {
		System.out.println("Evaluation loci: " + eval_num + "\tStandard loci: " + comp_num);
		System.out.println("Evaluation loci not called: " + eval_num_no_call + "\tEvaluation loci filtered out in false negative loci: " + eval_num_filtered);
		System.out.println("True Positive loci: " + tp);
		System.out.println("False Positive loci: " + fp);
		System.out.println("True Negative loci: " + tn);
		System.out.println("False Negative loci: " + fn);
		System.out.println("False Discoverary rate: " + (double)fp/(double)(fp+tp));
		System.out.println("Sensitivity rate: " + (double)tp/(double)(fn+tp));
	}
	
	private void comparison(VariantContext vc_comp,RefMetaDataTracker tracker){
		List<VariantContext> eval_bindings = tracker.getValues(eval);
		if(!eval_bindings.isEmpty()){
			VariantContext vc_eval = eval_bindings.get(0);
			
			eval_num++;
			if(homC){
				if(isHomoC(vc_comp) && vc_comp.getPhredScaledQual()>=mq){
					
					if(isHomoC(vc_eval) && vc_eval.getPhredScaledQual()>=mq){
						tp++;
					}
					else{
						fn++;
						if(vc_eval.getPhredScaledQual()<mq){
							eval_num_filtered++;
						}
					}
					
				}
				else{
					if(isHomoC(vc_eval) && vc_eval.getPhredScaledQual()>=mq){
						fp++;
					}
					else{
						tn++;
					}
				}
			}
			else{
				if(isHetSnp(vc_comp) && vc_comp.getPhredScaledQual()>=mq){
					if(isHetSnp(vc_eval) && vc_eval.getPhredScaledQual()>=mq){
						tp++;
					}
					else{
						fn++;
						if(vc_eval.getPhredScaledQual()<mq){
							eval_num_filtered++;
						}
					}
					
				}
				else{
					if(isHetSnp(vc_eval) && vc_eval.getPhredScaledQual()>=mq){
						fp++;
					}
					else{
						tn++;
					}
				}
			}
		}
		else{
			eval_num_no_call++;
			if(homC){
				if(isHomoC(vc_comp) && vc_comp.getPhredScaledQual()>=mq){
					fn++;
					
				}
				else{
					tn++;
				}
			}
			else{
				if(isHetSnp(vc_comp) && vc_comp.getPhredScaledQual()>=mq){
					fn++;
				}
				else{
					tn++;
				}
			}
		}
	}
	
	public boolean isHetSnp(VariantContext vc) {
     	 if(vc.hasGenotypes()){
     		Iterator<Genotype> it = vc.getGenotypes().iterator();
     		 while(it.hasNext()){
     			 if(it.next().isHet())
     				 return true;
     		 }
       }
          return false;
     }
	
	 public boolean isHomoC(VariantContext vc) {
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
}
