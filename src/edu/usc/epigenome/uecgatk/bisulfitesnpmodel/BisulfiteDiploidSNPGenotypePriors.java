package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.util.Arrays;
import java.util.Collection;
import java.util.EnumSet;

import org.broad.tribble.Feature;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.datasources.providers.RodLocusView;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypePriors;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;

public class BisulfiteDiploidSNPGenotypePriors implements GenotypePriors {
	// --------------------------------------------------------------------------------------------------------------
    //
    // Constants and static information
    //
    // --------------------------------------------------------------------------------------------------------------
    public static double HETEROZYGOSITY = 1e-3;
    
    public static double DBSNP_VALIDATE_HETEROZYGOSITY = 0.1;
    public static double DBSNP_NOVAL_HETEROZYGOSITY = 0.02;
    
    protected static double BISULFITE_CONVERSION_RATE = 0; 
    protected static double OVER_CONVERSION_RATE = 0;
    protected static double CYTOSINE_METHYLATION_RATE = 0;

    private final static double[] flatPriors = new double[DiploidGenotype.values().length];

    private double[] priors = null;

    /**
     * Create a new DiploidGenotypePriors object with flat priors for each diploid genotype
     */
    public BisulfiteDiploidSNPGenotypePriors() {
        priors = flatPriors.clone();
    }

    /**
     * Create a new GenotypeLikelihoods object with priors for a diploid with heterozygosity and reference
     * base ref
     *
     * @param ref
     * @param heterozygosity
     * @param probOfTriStateGenotype The prob of seeing a true B/C het when the reference is A
     */
    public BisulfiteDiploidSNPGenotypePriors(byte ref, double heterozygosity, double probOfTriStateGenotype) {
      
    	priors = flatPriors.clone();
    	HETEROZYGOSITY = heterozygosity;
    }
    
    public BisulfiteDiploidSNPGenotypePriors(byte ref, double heterozygosity, double probOfTriStateGenotype, double bisulfiteConversionRate) {
    	priors = flatPriors.clone();
    	HETEROZYGOSITY = heterozygosity;
        
    }

    public BisulfiteDiploidSNPGenotypePriors(double[] log10Priors) {
        priors = log10Priors.clone();
    }

    
    public double[] getPriors() {
        return priors;
    }
    
    //set prior for each locus
    public void setPriors(RefMetaDataTracker tracker, ReferenceContext ref, double heterozygosity, double probOfTriStateGenotype, double novelDbsnpHet, double validateDbsnpHet, GenomeLoc loc) {
    	
    	DBSNP_NOVAL_HETEROZYGOSITY = novelDbsnpHet;
    	DBSNP_VALIDATE_HETEROZYGOSITY = validateDbsnpHet;
    	
        byte refBase = ref.getBase();
        
        Collection<VariantContext> contexts = tracker.getVariantContexts(ref, DbSNPHelper.STANDARD_DBSNP_TRACK_NAME, null, loc, true, false);
        int count = 0;
        if(contexts != null){
        	for(VariantContext tmpVc : contexts){
                if(tmpVc.isSNP()){
                	priors = getReferencePolarizedPriors(refBase, DBSNP_VALIDATE_HETEROZYGOSITY, probOfTriStateGenotype);
                	count++;
                	break;
                }
                
            }
        }
        	if(count == 0)
        		priors = getReferencePolarizedPriors(refBase, heterozygosity, probOfTriStateGenotype);

    }
    
    
    /**
     * Returns the prior associated with DiploidGenotype g
     * @param g
     * @return log10 prior as a double
     */
    public double getPrior(DiploidGenotype g) {
        return getPriors()[g.ordinal()];
    }

    /**
     *
     * @param h
     * @return
     */
    public static double heterozygosity2HomRefProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        double v = 1.0 - (3 * h / 2.0);
        if (MathUtils.isNegative(v)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return v;
    }

    public static double heterozygosity2HetProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return h;
    }

    public static double heterozygosity2HomVarProbability(double h) {
        if (MathUtils.isNegative(h)) {
            throw new RuntimeException(String.format("Heterozygous value is bad %f", h));
        }

        return h / 2.0;
    }
    
    /**
    *check if the base is transition or transversion relative to reference base
    */
    public static boolean isTransition(byte base, byte ref) {
        boolean transition;
    	switch(ref){
        	case BaseUtils.A: transition = base == BaseUtils.G ? true : false; break;
        	case BaseUtils.G: transition = base == BaseUtils.A ? true : false; break;
        	case BaseUtils.C: transition = base == BaseUtils.T ? true : false; break;
        	case BaseUtils.T: transition = base == BaseUtils.C ? true : false; break;
        	default: throw new RuntimeException(String.format("ref base is bad "));
        }
    	return transition;
    }


    /**
     * Takes reference base, and three priors for hom-ref, het, hom-var, and fills in the priors vector
     * appropriately.
     *
     * Suppose A is the reference base, and we are given the probability of being hom-ref, het, and hom-var,
     * and that pTriSateGenotype (product of reference genome error and heterozygousity) is the true probability of observing reference A and a true genotype of B/C
     * then this sets the priors to:
     *
     * AA = pHomRef
     * AG = (pHet - pTriStateGenotype) * transitionRate
     * AC = AT = (pHet - pTriStateGenotype) * transversionRate
     * GG = pHomVar * transitionRate
     * CC = TT = pHomVar * transversionRate
     * CT = pTriStateGenotype * transversionRate
     * CG = GT = pTriStateGenotype * transitionRate
     * When i use switch clause, the program turns to be very slow.,so even this code is ugly for a lot of if, else, but it is much quicker...
     *
     * @param ref
     * @param heterozyosity
     * @param pRefError
     */

    public static double[] getReferencePolarizedPriors(byte ref, double heterozyosity, double pRefError) {
        if ( ! MathUtils.isBounded(pRefError, 0.0, 0.01) ) {
            throw new RuntimeException(String.format("BUG: p Reference error is out of bounds (0.0 - 0.01) is allow range %f", pRefError));
        }

        double pTriStateGenotype = heterozyosity * pRefError;

        double pHomRef = heterozygosity2HomRefProbability(heterozyosity);
        double pHet    = heterozygosity2HetProbability(heterozyosity);
        double pHomVar = heterozygosity2HomVarProbability(heterozyosity);

        if (MathUtils.compareDoubles(pHomRef + pHet + pHomVar, 1.0) != 0) {
            throw new RuntimeException(String.format("BUG: Prior probabilities don't sum to one => %f, %f, %f", pHomRef, pHet, pHomVar));
        }

        double[] priors = new double[DiploidGenotype.values().length];

        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            double POfG = 1;

            final double transitionRate = 2.0/3.0;
            final double transversionRate = 1.0/6.0;
           
            if(BaseUtils.basesAreEqual(ref, BaseUtils.A)){								//**ref = A
            	if(BaseUtils.basesAreEqual(g.base1, BaseUtils.A)){
            		if(BaseUtils.basesAreEqual(g.base2, BaseUtils.A)){
            			POfG = pHomRef;													//genotype: AA
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = (pHet - pTriStateGenotype) * transitionRate;				//genotype: AG
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.C)){
                    	POfG = (pHet - pTriStateGenotype) * transversionRate;			//genotype: AC
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate);			//genotype: AT
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.G)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = pHomVar * transitionRate;								//genotype: GG
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pTriStateGenotype * transitionRate;						//genotype: GT
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.C)){
                   if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = pTriStateGenotype * transitionRate;						//genotype: CG
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.C)){
                    	POfG = pHomVar * transversionRate;								//genotype: CC
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pTriStateGenotype * transversionRate;					//genotype: CT
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.T)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pHomVar * transversionRate;								//genotype: TT
                    }
                    else{
                    	
                    }
                }
                else{
                	
                }
            }
            else if(BaseUtils.basesAreEqual(ref, BaseUtils.G)){								//**ref = G
            	if(BaseUtils.basesAreEqual(g.base1, BaseUtils.A)){
            		if(BaseUtils.basesAreEqual(g.base2, BaseUtils.A)){
            			POfG = pHomVar * transitionRate;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = (pHet - pTriStateGenotype) * transitionRate;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.C)){
                    	POfG = (pTriStateGenotype * transitionRate);
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pTriStateGenotype * transitionRate;
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.G)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = pHomRef;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate);
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.C)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate);
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.C)){
                    	POfG = (pHomVar * transversionRate);
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = (pTriStateGenotype * transversionRate);
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.T)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pHomVar * transversionRate;
                    }
                    else{
                    	
                    }
                }
                else{
                	
                }
            }
            else if(BaseUtils.basesAreEqual(ref, BaseUtils.C)){								//**ref = C
            	if(BaseUtils.basesAreEqual(g.base1, BaseUtils.A)){
            		if(BaseUtils.basesAreEqual(g.base2, BaseUtils.A)){
            			POfG = pHomVar * transversionRate;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = pTriStateGenotype * transversionRate;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.C)){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate);
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pTriStateGenotype * transitionRate;
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.G)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = pHomVar * transversionRate;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pTriStateGenotype * transitionRate;
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.C)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = ((pHet - pTriStateGenotype) * transversionRate);
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.C)){
                    	POfG = pHomRef;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = ((pHet - pTriStateGenotype) * transitionRate);
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.T)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pHomVar * transitionRate;
                    }
                    else{
                    	
                    }
                }
                else{
                	
                }
            }
            else if(BaseUtils.basesAreEqual(ref, BaseUtils.T)){								//**ref = T
            	if(BaseUtils.basesAreEqual(g.base1, BaseUtils.A)){
            		if(BaseUtils.basesAreEqual(g.base2, BaseUtils.A)){
            			POfG = pHomVar * transversionRate;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = pTriStateGenotype * transversionRate;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.C)){
                    	POfG = (pTriStateGenotype * transitionRate);
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = (pHet - pTriStateGenotype) * transversionRate;
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.G)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.G)){
                    	POfG = pHomVar * transversionRate;
                    }
                    else if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = (pHet - pTriStateGenotype) * transversionRate;
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.C)){
                    if(BaseUtils.basesAreEqual(g.base2 , BaseUtils.G)){
                    	POfG = (pTriStateGenotype * transitionRate);
                    }
                    else if(BaseUtils.basesAreEqual(g.base2 , BaseUtils.C)){
                    	POfG = (pHomVar * transitionRate);
                    }
                    else if(BaseUtils.basesAreEqual(g.base2 , BaseUtils.T)){
                    	POfG = ((pHet - pTriStateGenotype) * transitionRate);
                    }
                    else{
                    	
                    }
                }
                else if(BaseUtils.basesAreEqual(g.base1, BaseUtils.T)){
                    if(BaseUtils.basesAreEqual(g.base2, BaseUtils.T)){
                    	POfG = pHomRef;
                    }
                    else{
                    	
                    }
                }
                else{
                	
                }
            }
            else{
            	
            }

            priors[g.ordinal()] = Math.log10(POfG);
        }

        return priors;
    }
  
    /*
     * made flat prior for each genotype
     */
    static {
        for ( DiploidGenotype g : DiploidGenotype.values() ) {
            flatPriors[g.ordinal()] = Math.log10(1.0 / DiploidGenotype.values().length);
        }
    }
	@Override
	public double getHeterozygosity() {
		// TODO Auto-generated method stub
		return HETEROZYGOSITY;
	}

	@Override
	public boolean validate(boolean throwException) {
		// TODO Auto-generated method stub
		return false;
	}

}
