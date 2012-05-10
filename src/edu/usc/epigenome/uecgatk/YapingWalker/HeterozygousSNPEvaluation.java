/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.util.Collection;
import java.util.Iterator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;


/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 5, 2012 9:00:51 PM
 * 
 */
@Analysis(description = "Assess site accuracy and sensitivity of calling heterozygous SNPs")
public class HeterozygousSNPEvaluation extends VariantEvaluator {

	// todo -- note this isn't strictly allele away.  It's really focused on sites.  A/T call at a validated A/G site is currently counted as a TP
    @DataPoint(description = "nComp") int nComp = 0;
    @DataPoint(description = "TP") int TP = 0;
    @DataPoint(description = "FP") int FP = 0;
    @DataPoint(description = "FN") int FN = 0;
    @DataPoint(description = "TN") int TN = 0;

    @DataPoint(description = "Sensitivity", format = "%.2f") double sensitivity = 0;
    @DataPoint(description = "Specificity", format = "%.2f") double specificity = 0;
    @DataPoint(description = "PPV", format = "%.2f") double PPV = 0;
    @DataPoint(description = "FDR", format = "%.2f") double FDR = 0;

    @DataPoint(description = "CompHetEvalNoCall") int CompHetEvalNoCall = 0;
    @DataPoint(description = "CompHetEvalFiltered") int CompHetEvalFiltered = 0;
    @DataPoint(description = "CompHetEvalHet") int CompHetEvalHet = 0;
    @DataPoint(description = "CompHetEvalHom") int CompHetEvalHom = 0;

    @DataPoint(description = "CompHomEvalNoCall") int CompHomEvalNoCall = 0;
    @DataPoint(description = "CompHomEvalFiltered") int CompHomEvalFiltered = 0;
    @DataPoint(description = "CompHomEvalHet") int CompHomEvalHet = 0;
    @DataPoint(description = "CompHomEvalHom") int CompHomEvalHom = 0;

    @DataPoint(description = "CompFiltered") int CompFiltered = 0;
    @DataPoint(description = "Eval and comp have different alleles") int nDifferentAlleleSites = 0;

    private static final boolean REQUIRE_IDENTICAL_ALLELES = false;

    private enum SiteStatus { NO_CALL, FILTERED, HET, HOM }

    // Counts of ValidationSiteStatus x CallSiteStatus
    final int[][] counts = new int[SiteStatus.values().length][SiteStatus.values().length];

    @Override public int getComparisonOrder() { return 2; }
    @Override public boolean enabled() { return true; }

    @Override
    public void finalizeEvaluation() {
        for ( SiteStatus x : SiteStatus.values() )
            CompFiltered += getCounts(SiteStatus.FILTERED, x);

        CompHetEvalNoCall = getCounts(SiteStatus.HET, SiteStatus.NO_CALL);
        CompHetEvalFiltered = getCounts(SiteStatus.HET, SiteStatus.FILTERED);
        CompHetEvalHet = getCounts(SiteStatus.HET, SiteStatus.HET);
        CompHetEvalHom = getCounts(SiteStatus.HET, SiteStatus.HOM);

        CompHomEvalNoCall = getCounts(SiteStatus.HOM, SiteStatus.NO_CALL);
        CompHomEvalFiltered = getCounts(SiteStatus.HOM, SiteStatus.FILTERED);
        CompHomEvalHet = getCounts(SiteStatus.HOM, SiteStatus.HET);
        CompHomEvalHom = getCounts(SiteStatus.HOM, SiteStatus.HOM);

        TP = CompHetEvalHet;
        FN = CompHetEvalNoCall + CompHetEvalFiltered + CompHetEvalHom;
        FP = CompHomEvalHet;
        TN = CompHomEvalNoCall + CompHomEvalFiltered + CompHomEvalHom;

        for ( SiteStatus x : SiteStatus.values() )
            for ( SiteStatus y : SiteStatus.values() )
                nComp += getCounts(x, y);

      //  if ( nComp != TP + FN + FP + TN + CompFiltered )
     //       throw new ReviewedStingException("BUG: nComp != TP + FN + FP + TN + CompFiltered!");

        sensitivity = (100.0 * TP) / (TP + FN);
        specificity = (TN+FP > 0) ? (100.0 * TN) / (TN + FP) : 100.0;
        PPV = (100.0 * TP) / (TP + FP);
        FDR = (100.0 * FP) / (FP + TP);
    }

    private int getCounts(SiteStatus comp, SiteStatus eval) {
        return counts[comp.ordinal()][eval.ordinal()];
    }

    @Override
    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( comp != null ) { // we only need to consider sites in comp
            if ( REQUIRE_IDENTICAL_ALLELES && (eval != null && haveDifferentAltAlleles(eval, comp)))
                nDifferentAlleleSites++;
            else {
                SiteStatus evalStatus = calcSiteStatus(eval);
                SiteStatus compStatus = calcSiteStatus(comp);
                counts[compStatus.ordinal()][evalStatus.ordinal()]++;
            }
        }

        return null; // we don't capture any interesting sites
    }

    //
    // helper routines
    //
    public SiteStatus calcSiteStatus(VariantContext vc) {
        if ( vc == null ) return SiteStatus.NO_CALL;
        if ( vc.isFiltered() ) return SiteStatus.FILTERED;
        if ( isHetSnp(vc) ) return SiteStatus.HET;
        if(vc.hasGenotypes()) return SiteStatus.HOM;  // 

        return SiteStatus.NO_CALL;//
    }



    public boolean haveDifferentAltAlleles(VariantContext eval, VariantContext comp) {
        Collection<Allele> evalAlts = eval.getAlternateAlleles();
        Collection<Allele> compAlts = comp.getAlternateAlleles();
        if ( evalAlts.size() != compAlts.size() ) {
            return true;
        } else {
            // same size => every alt from eval must be in comp
            for ( Allele a : evalAlts ) {
                if ( ! compAlts.contains(a) ) {
//                    System.out.printf("Different alleles: %s:%d eval=%s comp=%s\n\t\teval=%s\n\t\tcomp=%s%n",
//                            eval.getChr(), eval.getStart(), eval.getAlleles(), comp.getAlleles(), eval, comp);
                    return true;
                }
            }

            return false;
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
	

}