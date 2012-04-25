package edu.usc.epigenome.uecgatk.BisSNP;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextBuilder;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums.OUTPUT_MODE;

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

public class BisulfiteGenotyperEngine{

	private BisulfiteArgumentCollection BAC = null;
	
	private BisulfiteVariantCallContext bisulfiteVariantCallContext = null;
	
	private GenomeAnalysisEngine toolkit = null;
	
	//private UnifiedGenotyperEngine UGE = null;
	//cytosine pattern and their status
	//private ThreadLocal<HashMap<String,CytosineTypeStatus>> ctss = new ThreadLocal<HashMap<String,CytosineTypeStatus>>();
	//private ThreadLocal<CytosineTypeStatus> ctss = new ThreadLocal<CytosineTypeStatus>();
	
    // the model used for calculating genotypes
    private ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel> bglcms = new ThreadLocal<BisulfiteSNPGenotypeLikelihoodsCalculationModel>();

    // the allele frequency likelihoods (allocated once as an optimization)
//    private ThreadLocal<HashMap<String,double[]>> log10AlleleFrequencyPosteriors = new ThreadLocal<HashMap<String,double[]>>();

    // the priors object
 //   private GenotypePriors genotypePriors;

    // samples in input
  //  private Set<String> samples = new TreeSet<String>();

    // the various loggers and writers
 //   private Logger logger = null;

    // fasta reference reader to supplement the edges of the reference sequence for long reads
  //  private IndexedFastaSequenceFile referenceReader;

    // the standard filter to use for calls below the confidence threshold but above the emit threshold
    private static final Set<String> filter = new HashSet<String>(1);

 //   private static final int MISMATCH_WINDOW_SIZE = 20;
    
    private static boolean autoEstimateC = false;
    private static boolean secondIteration = false;

	protected double MAX_PHRED = 1000000;
	private int SOMATIC_STATS = 5;
	
	public static final String LOW_QUAL_FILTER_NAME = "LowQual";
	

	public BisulfiteGenotyperEngine(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, BisulfiteArgumentCollection BAC, GenomeAnalysisEngine toolkit, boolean autoEstimateC, boolean secondIteration) {
		this.BAC = BAC.clone();
		
		this.toolkit = toolkit;
		this.autoEstimateC = autoEstimateC;
		this.secondIteration = secondIteration;
		filter.add(LOW_QUAL_FILTER_NAME);
		//this.UGE = initializeUnifiedGenotypeEngine();
		// TODO Auto-generated constructor stub
		calculateLikelihoodsAndGenotypes(tracker, refContext, rawContext);
		
	}
	
	//public UnifiedGenotyperEngine initializeUnifiedGenotypeEngine(){
	//	return new UnifiedGenotyperEngine(toolkit, (UnifiedArgumentCollection)BAC);
	//}

	
	/*
	public void setCytosineTypeStatus(CytosineTypeStatus cts){
		HashMap<String,CytosineTypeStatus> cytosineStatuss = new HashMap<String,CytosineTypeStatus>();
		for(String sample : samples){
			cytosineStatuss.put(sample, cts.clone());
		}
		this.ctss.set(cytosineStatuss);
	}
	*/
	public BisulfiteVariantCallContext getBisulfiteVariantCallContext(){
		
		return bisulfiteVariantCallContext;
	}
	
	/**
     * Compute full BisulfiteVariantCallContext at a given locus.
     *
     * @param tracker    the meta data tracker
     * @param refContext the reference base
     * @param rawContext contextual information around the locus
     * @return the BisulfiteVariantCallContext object
     */
    public void calculateLikelihoodsAndGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext) {
        Map<String, AlignmentContext> stratifiedContexts = AlignmentContextUtils.splitContextBySampleName(rawContext);
        HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs = new HashMap<String, BisulfiteContextsGenotypeLikelihoods>();
        VariantContext vc = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.COMPLETE, BCGLs); //get likelihood and methylation pattern information from BisulfiteSNPGenotypeLikelihoodsCalculationModel
        if(rawContext.getBasePileup().getLocation().getStart() == BAC.testLocus){
        	
          	System.err.println("??" + vc.getAlleles().toString());
        
        }
        if ( vc == null )
            return;
        
      
        
        //VariantCallContext vcc = UGE.calculateGenotypes(tracker, refContext, rawContext, vc);
        bisulfiteVariantCallContext = calculateGenotypes(tracker, refContext, rawContext, stratifiedContexts, BCGLs, vc); //including all Reads group genotypes information into vcc, and provide most probable genotype and alt-allele for all of ReadsGroup 
 
        //return vcc;
    }

    

    public VariantContext calculateLikelihoods(RefMetaDataTracker tracker, ReferenceContext refContext, Map<String, AlignmentContext> stratifiedContexts, AlignmentContextUtils.ReadOrientation type, HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs) {
		if ( stratifiedContexts == null ){
			 return null;
		}

        // initialize the data for this thread if that hasn't been done yet
        if ( bglcms.get() == null ) {
            bglcms.set(BisulfiteGenotyperEngine.getGenotypeLikelihoodsCalculationObject(BAC, autoEstimateC, secondIteration));
           // if(refContext.getLocus().getStart() == BAC.testLocus){
            //	System.err.println("initiate\t" + refContext.getLocus().getStart());
           // }
        }

        BisulfiteSNPGenotypeLikelihoodsCalculationModel bglcm = (BisulfiteSNPGenotypeLikelihoodsCalculationModel) bglcms.get();
        //bglcm.initialize();
        
        bglcm.setBsLikelihoods(tracker, refContext, stratifiedContexts, type, BCGLs);
    
        
        if (!BCGLs.isEmpty()){
        	
        	return createVariantContextFromLikelihoods(refContext, bglcm.getRefAllele(), BCGLs);
        }
        else{
        	return null; 
        }
        
         
        
    }
    

    private void assignAFPosteriors(double[]likelihoods, double[] log10AFPosteriors){
    	for(int i = 0; i < likelihoods.length; i++){
    		log10AFPosteriors[i] = likelihoods[i];
    	}
    		
    }
	

    private BisulfiteVariantCallContext calculateGenotypes(RefMetaDataTracker tracker, ReferenceContext refContext, AlignmentContext rawContext, Map<String, AlignmentContext> stratifiedContexts, HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs, VariantContext vc) {
		// initialize the data for this thread if that hasn't been done yet
        if ( bglcms.get() == null ) {
            return null;
        }
       // System.err.println(vc.toString());
        HashMap<String, Object> attributes = new HashMap<String, Object>();
        int totalDepth = 0;
        int[] alleleCount = new int[2];
        alleleCount[1] = -1;
       Allele bestAlternateAlelle = Allele.NO_CALL;
        
        double logRatio = Double.NEGATIVE_INFINITY;
        int bestAF = 0;
        int numOfC = 0;
        int numOfT = 0;
        char Strand = '+';
        String cType = null;
        HashSet<String> cytosineConfirmed = new HashSet<String>();
        HashMap<String,Double> logRatiosMap = new HashMap<String,Double>();
        for ( String sample : BCGLs.keySet() ) {
        	BisulfiteContextsGenotypeLikelihoods GL = BCGLs.get(sample);
        	
        	totalDepth += GL.getDepth();
        	
        	 double[] log10AlleleFrequencyPosteriors = new double[3];
        	assignAFPosteriors(GL.getLikelihoods(),log10AlleleFrequencyPosteriors);
        	int bestAFguess = MathUtils.maxElementIndex(log10AlleleFrequencyPosteriors);
            int secondAFguess = MathUtils.minElementIndex(log10AlleleFrequencyPosteriors);
            for (int i = 0; i < log10AlleleFrequencyPosteriors.length; i++){
            	if(i != bestAFguess){	
            		if(log10AlleleFrequencyPosteriors[i] >= log10AlleleFrequencyPosteriors[secondAFguess]){
            			secondAFguess = i;
            		}
            	}
            }
           // double[] normalizedPosteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors, true, false);
            double[] normalizedPosteriors = log10AlleleFrequencyPosteriors;
            if(bestAFguess!=0){
            	bestAF = bestAFguess;
            	if(bestAlternateAlelle.isNoCall()){
            		bestAlternateAlelle = GL.getAlleleB();
            		alleleCount[0]+= bestAFguess==1?1:2;
            	}
            	else{
            		if(GL.getAlleleB().basesMatch(bestAlternateAlelle)){
            			alleleCount[0]+= bestAFguess==1?1:2;
            		}
            		else{
            			alleleCount[1]+= bestAFguess==1?1:2;
            		}
            	}
            	
            }
            	
            
            double logRatioTmp = 10 * (normalizedPosteriors[bestAFguess] - normalizedPosteriors[secondAFguess]);
            logRatiosMap.put(sample, logRatioTmp);
           // double logRatioTmp = 10 * (normalizedPosteriors[bestAFguess]-Math.log10(1- Math.pow(10, normalizedPosteriors[bestAFguess]) ));
            if(logRatioTmp > logRatio)
            	logRatio = logRatioTmp;
           // System.err.println(log10AlleleFrequencyPosteriors[bestAFguess] + "\t" + 10 * normalizedPosteriors[bestAFguess] + "\t" + logRatioTmp);
            if(passesCallThreshold(logRatioTmp) && GL.getBestMatchedCytosinePattern() != null){
            	numOfC += GL.getNumOfCReadsInBisulfiteCStrand();
            	numOfT += GL.getNumOfTReadsInBisulfiteCStrand();
            	Strand = GL.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).cytosineStrand;
            	for(String cytosinePattern : GL.getCytosineParameters().keySet()){
					
					if(GL.getCytosineParameters().get(cytosinePattern).isCytosinePattern){
						//if(cType == null){
						//	cType = cytosinePattern;
						//}
						//else{
						//	cType = cType + "," + cytosinePattern;
						//}
						if(GL.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).isHeterozygousCytosinePattern){
							cytosineConfirmed.add(BaseUtilsMore.makeIupacCodeFrom2String(GL.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).patternOfAlleleA, GL.getCytosineParameters().get(GL.getBestMatchedCytosinePattern()).patternOfAlleleB));
						}
						else{
							cytosineConfirmed.add(cytosinePattern);
						}
						
					}
            	}
            }
            
            if ( !BAC.NO_SLOD && bestAF != 0 ) {
                	//if(!passesCallThreshold(logRatioTmp)){ //only calculate SB score for Heterozygous SNP and Homozygous loci not pass threshold.. why?
                		// the overall lod
                        //double overallLog10PofNull = log10AlleleFrequencyPosteriors.get()[0];
                        double overallLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors, 1);
                        

                        // the forward lod
                        HashMap<String, BisulfiteContextsGenotypeLikelihoods> tmpGLs = new HashMap<String, BisulfiteContextsGenotypeLikelihoods>();
                        VariantContext vcForward = calculateLikelihoods(tracker, refContext, stratifiedContexts,AlignmentContextUtils.ReadOrientation.FORWARD, tmpGLs);
                        for ( int i = 0; i < log10AlleleFrequencyPosteriors.length; i++ )
                        	log10AlleleFrequencyPosteriors[i] = -1.0 * Double.MAX_VALUE;;
                        
                        if(tmpGLs.containsKey(sample)){
                        	assignAFPosteriors(tmpGLs.get(sample).getLikelihoods(),log10AlleleFrequencyPosteriors);
                        }
                        
                        //double[] normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
                        double forwardLog10PofNull = log10AlleleFrequencyPosteriors[0];
                        double forwardLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors, 1);
                       

                        // the reverse lod
                        tmpGLs = new HashMap<String, BisulfiteContextsGenotypeLikelihoods>();
                        VariantContext vcReverse = calculateLikelihoods(tracker, refContext, stratifiedContexts, AlignmentContextUtils.ReadOrientation.REVERSE, tmpGLs);
                        
                        for ( int i = 0; i < log10AlleleFrequencyPosteriors.length; i++ )
                        	log10AlleleFrequencyPosteriors[i] = -1.0 * Double.MAX_VALUE;;
                        
                        if(tmpGLs.containsKey(sample)){
                        	assignAFPosteriors(tmpGLs.get(sample).getLikelihoods(),log10AlleleFrequencyPosteriors);
                        }
                        
                        //normalizedLog10Posteriors = MathUtils.normalizeFromLog10(log10AlleleFrequencyPosteriors.get(), true);
                        double reverseLog10PofNull = log10AlleleFrequencyPosteriors[0];
                        double reverseLog10PofF = MathUtils.log10sumLog10(log10AlleleFrequencyPosteriors, 1);
                        

                        double forwardLod = forwardLog10PofF + reverseLog10PofNull - overallLog10PofF;
                        double reverseLod = reverseLog10PofF + forwardLog10PofNull - overallLog10PofF;
                        

                        // strand score is max bias between forward and reverse strands
                        double strandScore = Math.max(forwardLod, reverseLod);
                        
                      //  System.err.println(vcForward.getStart() + ": " +  forwardLog10PofF + "\t" + reverseLog10PofNull + "\t" + overallLog10PofF + "\t" + forwardLod + "\t" + reverseLog10PofF + "\t" + forwardLog10PofNull + "\t" + reverseLod + "\t" + strandScore);
                        // rescale by a factor of 10
                        strandScore *= 10.0;
                        //logger.debug(String.format("SLOD=%f", strandScore));

                        attributes.put("SB", Double.valueOf(strandScore));
                //	}
            }
        }  
        
            GenotypesContext genotypes = GenotypesContext.create();
            

            VariantContext dbsnp = null;

            // search for usable record
            //for( final VariantContext vc_input : tracker.getValues(BAC.dbsnp, refContext.getLocus()) ) {
            for( RODRecordList  rods: tracker.getBoundRodTracks() ) {
              //  if ( vc_input != null && vc_input.isSNP() ) {               
            	//if ( vc_input != null ) {    	
            	for(GATKFeature vc_input : rods){
            		if ( vc_input != null && vc_input.getUnderlyingObject() instanceof VariantContext) {
            			dbsnp = (VariantContext) vc_input.getUnderlyingObject();
            			break;
            		}
            	}
            }
            String rsID = null;
            if ( dbsnp != null && dbsnp.hasID()){
            	rsID = dbsnp.getID();
            	//System.err.println(dbsnp.getAttribute("VLD") + "\t" + rsID);
            //	 attributes.put(BisulfiteVCFConstants.ID_KEY, rsID);
            	 attributes.put(VCFConstants.DBSNP_KEY, true);
            }
           // System.err.println(BAC.dbsnp.toString() + "\t" + tracker.getBoundRodTracks().size());  

            // if the site was downsampled, record that fact
            if ( rawContext.hasPileupBeenDownsampled() )
                attributes.put(VCFConstants.DOWNSAMPLED_KEY, true);

            
            
            List<Allele> myAlleles = vc.getAlleles();
 
            GenomeLoc loc = refContext.getLocus();

            int endLoc = calculateEndPos(vc.getAlleles(), vc.getReference(), loc);
            //System.err.println(logRatio);
            
            
          //  int alleleNumber =  vc.getNAlleles();
          //  double[] alleleFrequency = new double[2];
            
            
            assignGenotypes(vc,bestAF, genotypes, logRatiosMap, BCGLs, BCGLs.keySet());
            
            attributes.put(BisulfiteVCFConstants.NUM_OF_SAMPLES, BCGLs.keySet().size());
            attributes.put(VCFConstants.DEPTH_KEY, totalDepth);
            /*
            if(bestAF != 0){
            	alleleFrequency[0] = (double)alleleCount[0]/(double)alleleNumber;
                if(alleleCount[1] == -1){
                	alleleFrequency[1] = (double)alleleCount[1]/(double)alleleNumber;
                    String alleleCountStr = String.format("%d",alleleCount[0]);
                    String alleleFreqStr = String.format("%.2f",alleleFrequency[0]);
                    attributes.put(VCFConstants.ALLELE_COUNT_KEY, alleleCountStr);
                    attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqStr);
                    attributes.put(VCFConstants.ALLELE_NUMBER_KEY, alleleNumber);
                }
                else{
                	alleleFrequency[1] = (double)alleleCount[1]/(double)alleleNumber;
                    String alleleCountStr = alleleCount[0] + "," + alleleCount[1];
                    String alleleFreqStr = String.format("%.2f",alleleFrequency[0]) + "," + String.format("%.2f", alleleFrequency[1]);
                    attributes.put(VCFConstants.ALLELE_COUNT_KEY, alleleCountStr);
                    attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, alleleFreqStr);
                    attributes.put(VCFConstants.ALLELE_NUMBER_KEY, alleleNumber);
                }
            	
            	
            }
            */
            
            //add INFO colum about methylation summary across Read Group
            if(passesCallThreshold(logRatio) & !cytosineConfirmed.isEmpty()){
            	for(String c : cytosineConfirmed){
            		if(cType == null){
						cType = c;
					}
					else{
						cType = cType + "," + c;
					}
            	}
            		//attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, numOfC);
            		//attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, numOfT);
                	attributes.put(BisulfiteVCFConstants.C_STRAND_KEY, Strand);
                	attributes.put(BisulfiteVCFConstants.CYTOSINE_TYPE, cType);
            }

            VariantContextBuilder vcb = new VariantContextBuilder("BG_call", loc.getContig(), loc.getStart(), endLoc,myAlleles);
            vcb.attributes(attributes);
            vcb.filters(passesCallThreshold(logRatio) ? null : filter);
            vcb.genotypes(genotypes);
            if ( rsID != null ){
            	vcb.id(rsID); 
            }
            	 
            vcb.log10PError(-logRatio/10.0);
            vcb.referenceBaseForIndel(refContext.getBase());
            VariantContext vcCall = vcb.make();
            BisulfiteVariantCallContext call = new BisulfiteVariantCallContext(BCGLs, vcCall, rawContext, refContext);
            call.confidentlyCalled = passesCallThreshold(logRatio);
            call.shouldEmit = passesEmitThreshold(logRatio);
            return call;
	}

    
	protected boolean passesEmitThreshold(double conf, int bestAFguess) {
        return (BAC.OutputMode != OUTPUT_MODE.EMIT_ALL_SITES || bestAFguess != 0) && conf >= Math.min(BAC.STANDARD_CONFIDENCE_FOR_CALLING, BAC.STANDARD_CONFIDENCE_FOR_EMITTING);
    }
	
		
	protected static BisulfiteSNPGenotypeLikelihoodsCalculationModel getGenotypeLikelihoodsCalculationObject(BisulfiteArgumentCollection BAC, boolean autoEstimateC, boolean secondIteration) {		
        	return new BisulfiteSNPGenotypeLikelihoodsCalculationModel(BAC, autoEstimateC, secondIteration);  	
    }

	
	protected boolean passesCallThreshold(double conf) {
        return conf >= BAC.STANDARD_CONFIDENCE_FOR_CALLING;
    }
	
	protected boolean passesEmitThreshold(double conf) {
        return conf >= BAC.STANDARD_CONFIDENCE_FOR_EMITTING;
    }
	
	private int calculateEndPos(List<Allele> alleles, Allele refAllele, GenomeLoc loc) {
        boolean isSNP = true;
        for (Allele a : alleles){
            if (a.getBaseString().length() != 1) {
                isSNP = false;
                break;
            }
        }

        int endLoc = loc.getStart();
        if ( !isSNP )
            endLoc += refAllele.length();

        return endLoc;
    }
	
	   /**
     * Can be overridden by concrete subclasses
     * @param vc                   variant context with genotype likelihoods
     * @param log10AlleleFrequencyPosteriors    allele frequency results
     * @param AFofMaxLikelihood    allele frequency of max likelihood
     *
     * @return calls
     */
	
	private void assignGenotypes(VariantContext vc,
                                                 int bestAF,
                                                 GenotypesContext calls,
                                                 HashMap<String,Double> logRatios,
                                                 HashMap<String,BisulfiteContextsGenotypeLikelihoods> BCGLs,
                                                 Set<String> sampleNames) {
        if ( !vc.isVariant() )
            throw new UserException("The VCF record passed in does not contain an ALT allele at " + vc.getChr() + ":" + vc.getStart());

        	GenotypesContext GLs = vc.getGenotypes();
        //	Set<String> sampleNames = GLs.getSampleNames();
        	for(String sampleName : sampleNames){
        		Genotype g = GLs.get(sampleName);
            	
             //   if ( !g.hasLikelihoods() )
             //       continue;
                Double logRatio = -logRatios.get(sampleName)/10.0;
                Allele alleleA = BCGLs.get(sampleName).getAlleleA();
                Allele alleleB = BCGLs.get(sampleName).getAlleleB();
                ArrayList<Allele> myAlleles = new ArrayList<Allele>();
                if (bestAF == 0) {
                    myAlleles.add(alleleA);
                    myAlleles.add(alleleA);
                    
                } else if(bestAF == 1) {
                    myAlleles.add(alleleA);
                    myAlleles.add(alleleB);
                   

                }  else {
                    myAlleles.add(alleleB);
                    myAlleles.add(alleleB);
                    
                }
               // System.err.println(sampleName + "\t" + myAlleles.toString() + "\t" + logRatio + "\t" + GLs.getSampleNames());
                //calls.put(sampleName, new Genotype(sampleName, myAlleles, logRatio, null, g.getAttributes(), false));
                calls.add(new Genotype(sampleName, myAlleles, logRatio, null, g.getAttributes(), false));
        	}
        	
        		
    }
	
    
    private VariantContext createVariantContextFromLikelihoods(ReferenceContext refContext, Allele refAllele, HashMap<String, BisulfiteContextsGenotypeLikelihoods> BCGLs) {
        
        List<Allele> noCall = new ArrayList<Allele>();
        noCall.add(Allele.NO_CALL);

        List<Allele> alleles = new ArrayList<Allele>();
        alleles.add(refAllele);
      //  boolean addedAltAllele = false;

        GenotypesContext genotypes = GenotypesContext.create();;
        for ( BisulfiteContextsGenotypeLikelihoods BCGL : BCGLs.values() ) {
         //   if ( !addedAltAllele ) {
         //       addedAltAllele = true;  
                if(!alleles.contains(BCGL.getAlleleA()))
                	alleles.add(BCGL.getAlleleA());
                if(!alleles.contains(BCGL.getAlleleB()))
                	alleles.add(BCGL.getAlleleB());
          //  }   
                HashMap<String, Object> attributes = new HashMap<String, Object>();
                GenotypeLikelihoods likelihoods = GenotypeLikelihoods.fromLog10Likelihoods(BCGL.getLikelihoods());
                attributes.put(VCFConstants.DEPTH_KEY, BCGL.getDepth());
                attributes.put(VCFConstants.GENOTYPE_POSTERIORS_KEY, likelihoods.getAsString());
                //attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, BCGL.getNumOfCReadsInBisulfiteCStrand());
                //attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, BCGL.getNumOfTReadsInBisulfiteCStrand()); //need to fix it to judge strand
                
                
                
                
                if(BCGL.getBestMatchedCytosinePattern() != null){
                	//attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, String.format("%.2f", 100*BCGL.getMethylationLevel()));
                	
                	attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, BCGL.getNumOfCReadsInBisulfiteCStrand());
            		attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, BCGL.getNumOfTReadsInBisulfiteCStrand());
                	if(BCGL.getCytosineParameters().get(BCGL.getBestMatchedCytosinePattern()).isHeterozygousCytosinePattern){
                		attributes.put(BisulfiteVCFConstants.BEST_C_PATTERN, BaseUtilsMore.makeIupacCodeFrom2String(BCGL.getCytosineParameters().get(BCGL.getBestMatchedCytosinePattern()).patternOfAlleleA, BCGL.getCytosineParameters().get(BCGL.getBestMatchedCytosinePattern()).patternOfAlleleB));
                	}
                	else{
                		attributes.put(BisulfiteVCFConstants.BEST_C_PATTERN, BCGL.getBestMatchedCytosinePattern());
                	}
            		attributes.put(BisulfiteVCFConstants.C_STATUS, BCGL.getBaseCountStatusAsString());
                }
                else{
                	attributes.put(BisulfiteVCFConstants.BEST_C_PATTERN, VCFConstants.MISSING_VALUE_v4);
                	//attributes.put(BisulfiteVCFConstants.CYTOSINE_METHY_VALUE, VCFConstants.MISSING_VALUE_v4);
                	attributes.put(BisulfiteVCFConstants.NUMBER_OF_C_KEY, VCFConstants.MISSING_VALUE_v4);
            		attributes.put(BisulfiteVCFConstants.NUMBER_OF_T_KEY, VCFConstants.MISSING_VALUE_v4);
                	attributes.put(BisulfiteVCFConstants.C_STATUS, VCFConstants.MISSING_VALUE_v4);
                }
                attributes.put(BisulfiteVCFConstants.SOMATIC_STAT_VAR, SOMATIC_STATS);
                attributes.put(BisulfiteVCFConstants.READS_SUPPORT_ALT, BCGL.getDP4AsString());
                attributes.put(BisulfiteVCFConstants.AVE_BASE_QUALITY_KEY, BCGL.getAveBaseQualAsString()); // need to ask kothiyalp@mail.nih.gov the average mean RMS or arithmic average here
                //Reads supporting ALT. Number of 1) forward ref alleles; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref alleles"
                //Average base quality for reads supporting alleles. For each allele, in the same order as listed
                genotypes.add(new Genotype(BCGL.getSample(), alleles, Genotype.NO_LOG10_PERROR, null, attributes, false));
                
        }

        GenomeLoc loc = refContext.getLocus();
        int endLoc = calculateEndPos(alleles, refAllele, loc);
       // System.err.println(alleles.toString());
        VariantContextBuilder vcb = new VariantContextBuilder("BG_call", loc.getContig(), loc.getStart(), endLoc,alleles);
        vcb.referenceBaseForIndel(refContext.getBase());
   	 	vcb.log10PError(VariantContext.NO_LOG10_PERROR);
   	 	vcb.genotypes(genotypes);
   	 	if(loc.getStart() == BAC.testLocus)
   	 		System.err.println("hah" + alleles.toString());
        return vcb.make();
    }
    
    
    
    
}
