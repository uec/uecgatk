/**
 * This is to implement Robert Shoemaker 2010 Genome Research paper's algorithm's comparison with BisSNP
 * it is only dealing with ShoemakerReadsConversionWalker's resul. It only looks at dbSNP position with 10X, then it will do overlapped with samtools result..
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Gather;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.filters.BadMateFilter;
import org.broadinstitute.sting.gatk.filters.NotPrimaryAlignmentFilter;
import org.broadinstitute.sting.gatk.filters.UnmappedReadFilter;
import org.broadinstitute.sting.gatk.io.StingSAMFileWriter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;
import org.broadinstitute.sting.gatk.walkers.recalibration.CountCovariatesGatherer;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.genomeLibs.FisherExactTest;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteGenotyper.ContextCondition;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObject;
import edu.usc.epigenome.uecgatk.YapingWriter.bedObjectWriterImp;
import jsc.contingencytables.FishersExactTest;
import jsc.contingencytables.ContingencyTable2x2;
/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time Apr 29, 2012 3:43:33 PM
 * 
 */
@WalkerName("shoemakerGenotyper")
@Reference(window=@Window(start=-200,stop=200))
@By(DataSource.REFERENCE)
@ReadFilters( {UnmappedReadFilter.class, BadMateFilter.class, NotPrimaryAlignmentFilter.class} )
public class shoemakerGenotyper extends LocusWalker<Long,Long> implements
		TreeReducible<Long> {

	 @Input(fullName="dbsnp", shortName = "D", doc="dbSNP file", required=false)
	 public RodBinding<VariantContext> dbsnp;
	 
	 @Output(fullName="output", shortName = "o", doc = "The output of genotype result", required = true)
	 public String file=null;
	 
	 @Argument(fullName = "qual_threshold", shortName = "qual",doc = "The genotype qauality score threshold", required = false)
	 public double qual =10;
	 
	 @Argument(fullName = "mimimum_coverage", shortName = "cov",doc = "Minimum coverage required", required = false)
	 public int cov =10;
	 

	 private int MIN_BASE_QUAL = 15;
	 private int MIN_COVERAGE = cov;
	 private double MIN_PERCENT_READS_REVERSE_STRAND = 0.2;
	 private int MAX_COV=250;
	 
	 private double MAX_INFI=-100000.0;
	 
	 private bedObjectWriterImp writer=null;
	 
	 private int deugPos=-1;
	 private static boolean verbose=false;
	 
	/**
	 * 
	 */
	public void initialize() {
		writer=new bedObjectWriterImp(new File(file));
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public Long map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		
		if(!context.hasBasePileup()){
			return null;
		}
		if(context.getBasePileup().depthOfCoverage() < MIN_COVERAGE || context.getBasePileup().depthOfCoverage() > MAX_COV){
			return null;
		}
		if(context.getBasePileup().getLocation().getStart()==deugPos){
			verbose=true;
		}
		else{
			verbose=false;
		}
		
		
		//System.err.println(context.getBasePileup().getLocation().getStart());
		ReadBackedPileup filteredPileup = getFilteredPileup(context.getBasePileup());
		if(filteredPileup == null || filteredPileup.depthOfCoverage()==0)
			return null;
		
		int[] nucleotideFrequencyMatrix = generateNucleotideMatrix(filteredPileup); // get normalized A,G,C,T which involved base quality information to attenuate bad base effect. 
		if(filteredPileup.getLocation().getStart()==deugPos){
			for(int value : nucleotideFrequencyMatrix){
				System.err.println(value);
			}
			
		}
			
		double[] genotypeProbabilityMatrix = getPvalueFromNucleotideFrequencyMatrix(nucleotideFrequencyMatrix, filteredPileup);
		DiploidGenotype bestGenotypeIndex = DiploidGenotype.createHomGenotype(BaseUtils.A);;
		DiploidGenotype secondBestGenotypeIndex = DiploidGenotype.createHomGenotype(BaseUtils.A);;
		double maximum= Double.NEGATIVE_INFINITY;
		double second= Double.NEGATIVE_INFINITY;
		for ( DiploidGenotype g : DiploidGenotype.values() ){
			if(filteredPileup.getLocation().getStart()==deugPos)
				System.err.println(g.toString() + "\t" + genotypeProbabilityMatrix[g.ordinal()]);
			if(genotypeProbabilityMatrix[g.ordinal()]>maximum){
				second=maximum;
				maximum=genotypeProbabilityMatrix[g.ordinal()];
				secondBestGenotypeIndex=bestGenotypeIndex;
				bestGenotypeIndex=g;
			}
			else if(genotypeProbabilityMatrix[g.ordinal()]>second){
				second=genotypeProbabilityMatrix[g.ordinal()];
				secondBestGenotypeIndex=g;
			}
		}
		double score = genotypeProbabilityMatrix[bestGenotypeIndex.ordinal()]-genotypeProbabilityMatrix[secondBestGenotypeIndex.ordinal()];
		
		if( 10*score > qual){
			List<Object> values = new ArrayList<Object>(); 
			values.add(bestGenotypeIndex.toString());
			values.add(secondBestGenotypeIndex.toString());
			values.add(10*score);
			
			List<VariantContext> comp_bindings = tracker.getValues(dbsnp);
			if(!comp_bindings.isEmpty()){
				VariantContext vc_comp = comp_bindings.get(0); 
				if ( vc_comp != null && vc_comp.hasID()){
					values.add(vc_comp.getID());
				}
				else{
					values.add(".");
				}
			}
			else{
				values.add(".");
			}
			
			bedObject obj = new bedObject(filteredPileup.getLocation().getContig(), filteredPileup.getLocation().getStart()-1, filteredPileup.getLocation().getStop(), values);
			writer.add(obj);
			//genotypeFile.printf("%s\t%d\t%s\t%s\t%.2f\n", filteredPileup.getLocation().getContig(), filteredPileup.getLocation().getStart(), bestGenotypeIndex.toString(), secondBestGenotypeIndex.toString(), score);  //chr, loc, bestGenotype,secondBestGenotype, qual
		}
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public Long reduceInit() {
		// TODO Auto-generated method stub
		return null;
	}



	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Long treeReduce(Long lhs, Long rhs) {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public Long reduce(Long value, Long sum) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void onTraversalDone(Long sum) {
		writer.close();
		logger.info("Finished!..");
	}

	
	/*
	 * check base and +/- 3 bases, all of them's quality score should be more than 15, return pileup only contain filtered bases which satisfy these criteria
	 * check loci, when good base > 20 % in one strand, its complement base should be > 20% in the other strand
	 */
	private ReadBackedPileup getFilteredPileup(ReadBackedPileup pileup){
		List<GATKSAMRecord> reads =  new ArrayList<GATKSAMRecord>();;
		List<Integer> elementOffsets = new ArrayList<Integer>();

		int[] baseFwd = new int[5]; // record A,C,G,T, total count(fwd) in fwd strand
		int[] baseRev = new int[5];// record A,C,G,T, total count(rev) in reverse strand
		int baseCountToCheckUpDownStream = 3;
		for ( PileupElement p : pileup ) {
			
			int i = 0-baseCountToCheckUpDownStream;
			byte[] quals = p.getRead().getBaseQualities();
			boolean goodAdjacentBase=true;
			while(i <= baseCountToCheckUpDownStream){
				int elementOffset = i + p.getOffset();
				if(elementOffset < 0 || elementOffset > p.getRead().getReadLength()-1){
					i++;
					continue;
				}
				if(quals[elementOffset] < MIN_BASE_QUAL){
					i++;
					goodAdjacentBase=false;
					continue;
				}
				if(pileup.getLocation().getStart()==deugPos){
					
				//	System.err.println(i + "\t" + elementOffset + "\t" + p.getOffset());
				
				}
				if(i==0){
					if(p.getRead().getReadNegativeStrandFlag()){
						baseRev[4]++;
						if(BaseUtils.basesAreEqual(p.getBase(), BaseUtils.A)){
							baseRev[0]++;
						}
						else if(BaseUtils.basesAreEqual(p.getBase(), BaseUtils.C)){
							baseRev[1]++;
						}
						else if(BaseUtils.basesAreEqual(p.getBase(), BaseUtils.G)){
							baseRev[2]++;
						}
						else if(BaseUtils.basesAreEqual(p.getBase(), BaseUtils.T)){
							baseRev[3]++;
						}
					}
					else{
						baseFwd[4]++;
						if(BaseUtils.basesAreEqual(p.getBase(), BaseUtils.A)){
							baseFwd[0]++;
						}
						else if(BaseUtils.basesAreEqual(p.getBase(), BaseUtils.C)){
							baseFwd[1]++;
						}
						else if(BaseUtils.basesAreEqual(p.getBase(), BaseUtils.G)){
							baseFwd[2]++;
						}
						else if(BaseUtils.basesAreEqual(p.getBase(), BaseUtils.T)){
							baseFwd[3]++;
						}
					}
				}
				
				i++;
			}
			if(!goodAdjacentBase){
				continue;
			}
			
			elementOffsets.add(p.getOffset());
			reads.add(p.getRead());
		}
		
		for(int i = 0; i< baseFwd.length-1;i++){
			if(pileup.getLocation().getStart()==deugPos){
				
				System.err.println(baseFwd[i] + "\t" + baseFwd[4] + "\tfwd:" + i + "\t" + baseRev[i] + "\t" + baseRev[4]);
			
			}
			if((double)baseFwd[i]/(double)baseFwd[4] >MIN_PERCENT_READS_REVERSE_STRAND){
				if((double)baseRev[i]/(double)baseRev[4] < MIN_PERCENT_READS_REVERSE_STRAND){
					return null;
				}
			}
		}
		for(int i = 0; i< baseRev.length-1;i++){
			if(pileup.getLocation().getStart()==deugPos){
				
				System.err.println(baseRev[i] + "\t" + baseRev[4] + "\trew:" + i);
			
			}
			if((double)baseRev[i]/(double)baseRev[4] >MIN_PERCENT_READS_REVERSE_STRAND){
				if((double)baseFwd[i]/(double)baseFwd[4] < MIN_PERCENT_READS_REVERSE_STRAND){
					return null;
				}
			}
		}
		
		return new ReadBackedPileupImpl(pileup.getLocation(),reads,elementOffsets);
	}
	
	private int[] generateNucleotideMatrix(ReadBackedPileup pileup){
		double[] nucleotideProbabilityMatrix = new double[4]; //A,C,G,T,uncertain
		for(int i=0; i<nucleotideProbabilityMatrix.length; i++){
			nucleotideProbabilityMatrix[i]=Double.NEGATIVE_INFINITY;
		}

		for ( PileupElement p : pileup ) {
			add(nucleotideProbabilityMatrix, p);
		}
		
		double[] nucleotideProbabilityMatrixNormalized = normalizeToLog10(nucleotideProbabilityMatrix, false, false, true);
		int[] nucleotideFrequencyMatrix = new int[4];
		for ( int i= 0; i < nucleotideProbabilityMatrixNormalized.length; i++ ) {
			if(Double.isNaN(nucleotideProbabilityMatrixNormalized[i])){
				nucleotideFrequencyMatrix[i] = 0;
			}
			else{
				nucleotideFrequencyMatrix[i] = (int) (nucleotideProbabilityMatrixNormalized[i] * pileup.depthOfCoverage());  
			}
			if(pileup.getLocation().getStart()==deugPos){
				System.err.println(nucleotideFrequencyMatrix[i] + "\t" + nucleotideProbabilityMatrixNormalized[i] + "\t" + pileup.depthOfCoverage() + "\t" + nucleotideProbabilityMatrix[i]);
			}
			
		}
		return nucleotideFrequencyMatrix;
	}
	
	private void add(double[] nucleotideProbabilityMatrix, PileupElement p){
		double probability = 1.0-Math.pow(10.0,(-p.getQual()/10.0));
		byte base = p.getBase();
		
		if(p.getRead().getReadNegativeStrandFlag()){
			//base=BaseUtils.simpleComplement(base);
		}
			if(BaseUtils.basesAreEqual(base, BaseUtils.A) ){   
		
				if(Double.isInfinite(nucleotideProbabilityMatrix[0])){
					nucleotideProbabilityMatrix[0] = probability;
				}
				else{
					nucleotideProbabilityMatrix[0] += probability;
				}
				
			}
			else if(BaseUtils.basesAreEqual(base, BaseUtils.C) ){  
				if(Double.isInfinite(nucleotideProbabilityMatrix[1])){
					nucleotideProbabilityMatrix[1] = probability;
				}
				else{
					nucleotideProbabilityMatrix[1] += probability;
				}
				//nucleotideProbabilityMatrix[1] += Math.log10(probability);
			}
			else if(BaseUtils.basesAreEqual(base, BaseUtils.G) ){ 
				if(Double.isInfinite(nucleotideProbabilityMatrix[2])){
					nucleotideProbabilityMatrix[2] = probability;
				}
				else{
					nucleotideProbabilityMatrix[2] += probability;
				}
				//nucleotideProbabilityMatrix[2] += Math.log10(probability);
			}
			else if(BaseUtils.basesAreEqual(base, BaseUtils.T) ){ 
				if(Double.isInfinite(nucleotideProbabilityMatrix[3])){
					nucleotideProbabilityMatrix[3] = probability;
				}
				else{
					nucleotideProbabilityMatrix[3] += probability;
				}
				//nucleotideProbabilityMatrix[3] += Math.log10(probability);
			}
			//else{
			//	nucleotideProbabilityMatrix[4] += Math.log10(probability);
			//}
			
			if(verbose){
				for(int i = 0; i<nucleotideProbabilityMatrix.length; i++){
					System.err.println(nucleotideProbabilityMatrix[i]);
				}
				System.err.println(probability);
			}
				
	}
	
	//return normalize log 10 p vlaue for each genotypes
	private double[] getPvalueFromNucleotideFrequencyMatrix(int[] nucleotideFrequencyMatrix, ReadBackedPileup pileup){
		double[] genotypeProbability = new double[DiploidGenotype.values().length];
		int coverage = pileup.depthOfCoverage();
		for ( DiploidGenotype g : DiploidGenotype.values() ) {
			int base1Count=getFrequencyFromMatrix(g.base1, nucleotideFrequencyMatrix);
			int base2Count=getFrequencyFromMatrix(g.base2, nucleotideFrequencyMatrix); 
			int base1OtherCount=getOtherBaseFrequencyFromMatrix(g.base1, nucleotideFrequencyMatrix);
			int base2OtherCount=getOtherBaseFrequencyFromMatrix(g.base2, nucleotideFrequencyMatrix);
		//	System.err.println(pileup.getLocation().getStart());
			if(g.isHet()){
				genotypeProbability[g.ordinal()]=getFisherPvalue(base1Count,base1OtherCount, (int)(coverage*0.5), (int)(coverage*0.5));
				genotypeProbability[g.ordinal()]+=getFisherPvalue(base2Count,base2OtherCount, (int)(coverage*0.5), (int)(coverage*0.5));
				
			}
			else{
				genotypeProbability[g.ordinal()]=getFisherPvalue(base1Count,base1OtherCount, coverage, 0);
				genotypeProbability[g.ordinal()]+=getFisherPvalue(base2Count,base2OtherCount, coverage, 0);
			}
			if(pileup.getLocation().getStart()==deugPos){
				
				System.err.println(pileup.getLocation().getStart() + "\t" + g.toString() + "\t" + getFrequencyFromMatrix(g.base1, nucleotideFrequencyMatrix) + "\t" + getOtherBaseFrequencyFromMatrix(g.base1, nucleotideFrequencyMatrix) + "\t" + getFrequencyFromMatrix(g.base2, nucleotideFrequencyMatrix) + "\t" + getOtherBaseFrequencyFromMatrix(g.base2, nucleotideFrequencyMatrix));
			
		}
		}
		return MathUtils.normalizeFromLog10(genotypeProbability, true, false);
	}
	
	private int getFrequencyFromMatrix(byte base, int[] nucleotideFrequencyMatrix){
		if(BaseUtils.basesAreEqual(base, BaseUtils.A) ){
			return nucleotideFrequencyMatrix[0];
		}
		else if(BaseUtils.basesAreEqual(base, BaseUtils.C)){
			return nucleotideFrequencyMatrix[1];
		}
		else if(BaseUtils.basesAreEqual(base, BaseUtils.G)){
			return nucleotideFrequencyMatrix[2];
		}
		else if(BaseUtils.basesAreEqual(base, BaseUtils.T)){
			return nucleotideFrequencyMatrix[3];
		}
		return -1;
	//	else{
	//		return nucleotideFrequencyMatrix[4];
	//	}
	}
	
	private int getOtherBaseFrequencyFromMatrix(byte base, int[] nucleotideFrequencyMatrix){
		if(BaseUtils.basesAreEqual(base, BaseUtils.A) ){
			return nucleotideFrequencyMatrix[1] + nucleotideFrequencyMatrix[2] + nucleotideFrequencyMatrix[3];
		}
		else if(BaseUtils.basesAreEqual(base, BaseUtils.C)){
			return nucleotideFrequencyMatrix[0] + nucleotideFrequencyMatrix[2] + nucleotideFrequencyMatrix[3];
		}
		else if(BaseUtils.basesAreEqual(base, BaseUtils.G)){
			return nucleotideFrequencyMatrix[0] + nucleotideFrequencyMatrix[1] + nucleotideFrequencyMatrix[3];
		}
		else if(BaseUtils.basesAreEqual(base, BaseUtils.T)){
			return nucleotideFrequencyMatrix[0] + nucleotideFrequencyMatrix[1] + nucleotideFrequencyMatrix[2];
		}
		else{
			return -1;
		//	return nucleotideFrequencyMatrix[0] + nucleotideFrequencyMatrix[1] + nucleotideFrequencyMatrix[2] + nucleotideFrequencyMatrix[3];
		}
	}
	
	//return log10 p value
	private double getFisherPvalue(int baseCount, int otherBaseCount, int expectBaseCount, int expectOtherBaseCount){
	//	System.err.println(baseCount + "\t" + otherBaseCount + "\t" + expectBaseCount + "\t" + expectOtherBaseCount);
	//	ContingencyTable2x2 table = new ContingencyTable2x2(baseCount, otherBaseCount, expectBaseCount, expectOtherBaseCount);
	//	ContingencyTable2x2 table = new ContingencyTable2x2(22, 0, expectBaseCount, expectOtherBaseCount);
	//	FishersExactTest fisherExact = new FishersExactTest(table);
	//	return Math.log10(fisherExact.getSP());
		FisherExactTest fisherExact = new FisherExactTest(1000);
		 double pValue = Math.log10(fisherExact.getTwoTailedP(baseCount,otherBaseCount,expectBaseCount,expectOtherBaseCount));
		// System.err.println(22 + "\t" + 0 + "\t" + expectBaseCount + "\t" + expectOtherBaseCount + "\t" + pValue);
		  return pValue;
	}
	
	public static double[] normalizeToLog10(double[] array, boolean takeLog10OfOutput, boolean keepInLogSpace, boolean ignoreInfiPos) {

        // for precision purposes, we need to add (or really subtract, since they're
        // all negative) the largest value; also, we need to convert to normal-space.
        double maxValue = Utils.findMaxEntry(array);
      //  if(verbose)
		//	System.err.println(maxValue);
        // we may decide to just normalize in log space without converting to linear space
        

        // default case: go to linear space
        double[] normalized = new double[array.length];

        for (int i = 0; i < array.length; i++){
        	normalized[i] = Math.pow(10,Math.log10(array[i]) - Math.log10(maxValue));
        	//if(verbose)
    		//	System.err.println("ha" + "\t" + normalized[i] + "\t" + array[i] + "\t" + maxValue);
        }
            

        // normalize
        double sum = 0.0;
        for (int i = 0; i < array.length; i++){
        	if(!Double.isNaN(normalized[i])){
        		sum += normalized[i];
        	}
        }
            
        for (int i = 0; i < array.length; i++) {
            double x = normalized[i] / sum;
            if (takeLog10OfOutput)
                x = Math.log10(x);
            normalized[i] = x;
        }

        return normalized;
    }

}
