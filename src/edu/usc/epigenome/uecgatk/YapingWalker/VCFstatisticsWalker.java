/**
 * 
 */
package edu.usc.epigenome.uecgatk.YapingWalker;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

import net.sf.samtools.SAMSequenceDictionary;

import org.apache.commons.math.stat.descriptive.summary.Sum;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.BisSNP.BisSNPUtils;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteGenotyper;
import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteVCFConstants;
import edu.usc.epigenome.uecgatk.BisSNP.SortingTcgaVCFWriter;
import edu.usc.epigenome.uecgatk.BisSNP.TcgaVCFWriter;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time May 11, 2012 12:26:50 PM
 * 
 */
@Reference(window=@Window(start=-200,stop=200))
public class VCFstatisticsWalker extends RodWalker<VCFstatisticsWalker.VCFplusRef, VCFstatisticsWalker.VCFCondition>  implements TreeReducible<VCFstatisticsWalker.VCFCondition>  {

	/**
     * The variant file(s) to evaluate.
     */
    @Input(fullName="eval", shortName = "eval", doc="Input evaluation file(s)", required=true)
    public RodBinding<VariantContext> eval;
    
    @Output(doc="Output summary statistics", required=true)
    public String outFile;
    
    @Output(fullName="RefCph", shortName = "RefCph", doc="Output summary statistics", required=false)
    public String outRefCphFile;
    
    @Output(fullName="RefCpg", shortName = "RefCpg", doc="Output summary statistics", required=false)
    public String outRefCpgFile;
    
    @Argument(fullName="genotype_qual", shortName = "qual", doc="genotype quality score filter for heterozygous SNP, default: 20", required=false)
    public double qual=20;
    
    @Argument(fullName="strand_bias", shortName = "sb", doc="strand bias filter for heterozygous SNP, default: 0", required=false)
    public double sb=0;
    
    @Argument(fullName="max_coverage", shortName = "maxCov", doc="maximum coverage filter for heterozygous SNP, default: 100000", required=false)
    public int maxCov=100000;
    
    @Argument(fullName="quality_by_depth", shortName = "qd", doc="quality by depth filter for heterozygous SNP, default: 0", required=false)
    public double qd=0;
    
    @Argument(fullName="mapping_quality_zero", shortName = "mq0", doc="fraction of mapping_quality_zero filter for heterozygous SNP, default: 1.0", required=false)
    public double mq0=1.0;
    
    @Argument(fullName="min_ct_coverage", shortName = "minCT", doc="minimum number of CT reads for count methylation level, default: 1", required=false)
    public int minCT=1;
    
    @Argument(fullName="min_bq", shortName = "minBQ", doc="minimum base quality for both of strand, default: 10, not use this option yet", required=false)
    public int minBQ=10;
	
    private PrintStream writer = null;
    protected TcgaVCFWriter refCphWriter = null;
    protected TcgaVCFWriter refCpgWriter = null;
    
    protected SortingTcgaVCFWriter multiThreadRefCphWriter = null;
    protected SortingTcgaVCFWriter multiThreadRefCpgWriter = null;
    
    private int MAXIMUM_CACHE_FOR_OUTPUT_VCF = 10000000;
    
    public static class VCFCondition {
    	long nLociVisited = 0;
    	
    	long nLociFilterOut = 0;
    	//cytosine number
    	long nCytosineTotal = 0;
    	
    //	double sumCytosineTotal = 0;
    	long nCpgTotal = 0;
    	
    	//long[][][] num = new long[2][2][3];   //Cpg(0) or Cph(1); Hom(0) or Het(1); RefCpg(0) or RefCph(1) or RefNotC(2)
    			
    //	double[][][] methy =  new double[2][2][3]; 
    	
    	long nCpgHom = 0;
    	
    	long nCpgHomRefCpg = 0;
    	
    	long nCpgHomRefCph = 0;
    	
    	long nCpgHomRefNotC = 0;
    	
    	long nCpgHet = 0;
    	
    	long nCpgHetRefCpg = 0;
    	
    	long nCpgHetRefCph = 0;
    	
    	long nCpgHetRefNotC = 0;
    	
    	long nCpgHetAtC_RefCpg = 0;
    	
    	long nCpgHetAtC_RefCph = 0;
    	
    	long nCpgHetAtC_RefNotC = 0;
    	
    	long nCpgHetAtG_RefCpg = 0;
    	
    	long nCpgHetAtG_RefCph = 0;
    	
    	long nCpgHetAtG_RefNotC = 0;
    	
    	long nCphTotal = 0;

    	long nCphHom = 0;
    	
    	long nCphHomRefCpg = 0;
    	
    	long nCphHomRefCph = 0;
    	
    	long nCphHomRefNotC = 0;
    	
    	long nCphHet = 0;
    	
    	long nCphHetRefCpg = 0;
    	
    	long nCphHetRefCph = 0;
    	
    	long nCphHetRefNotC = 0;
    	
    	long nCphHetAtC_RefCpg = 0;
    	
    	long nCphHetAtC_RefCph = 0;
    	
    	long nCphHetAtC_RefNotC = 0;
    	

    	
    	//methy level
    	
    	double sumCytosineTotal = 0;
    	
    	double sumCpgTotal = 0;
    	
    	double sumCpgHom = 0;
    	
    	double sumCpgHomRefCpg = 0;
    	
    	double sumCpgHomRefCph = 0;
    	
    	double sumCpgHomRefNotC = 0;
    	
    	double sumCpgHet = 0;
    	
    	double sumCpgHetRefCpg = 0;
    	
    	double sumCpgHetRefCph = 0;
    	
    	double sumCpgHetRefNotC = 0;
    	
    	double sumCpgHetAtC_RefCpg = 0;
    	
    	double sumCpgHetAtC_RefCph = 0;
    	
    	double sumCpgHetAtC_RefNotC = 0;
    	
    	double sumCpgHetAtG_RefCpg = 0;
    	
    	double sumCpgHetAtG_RefCph = 0;
    	
    	double sumCpgHetAtG_RefNotC = 0;
    	
    	double sumCphTotal = 0;

    	double sumCphHom = 0;
    	
    	double sumCphHomRefCpg = 0;
    	
    	double sumCphHomRefCph = 0;
    	
    	double sumCphHomRefNotC = 0;
    	
    	double sumCphHet = 0;
    	
    	double sumCphHetRefCpg = 0;
    	
    	double sumCphHetRefCph = 0;
    	
    	double sumCphHetRefNotC = 0;
    	
    	double sumCphHetAtC_RefCpg = 0;
    	
    	double sumCphHetAtC_RefCph = 0;
    	
    	double sumCphHetAtC_RefNotC = 0;
    	
    }
    
	
	/**
	 * 
	 */
	public void initialize() {
		try {
			writer = new PrintStream(new File(outFile));
			if(outRefCphFile != null && outRefCpgFile != null){
				SAMSequenceDictionary refDict = getToolkit().getMasterSequenceDictionary();
				refCphWriter = new TcgaVCFWriter(new File(outRefCphFile), refDict);
				refCpgWriter = new TcgaVCFWriter(new File(outRefCpgFile), refDict);
				if(getToolkit().getArguments().numberOfThreads>1){
					multiThreadRefCphWriter = new SortingTcgaVCFWriter(refCphWriter, MAXIMUM_CACHE_FOR_OUTPUT_VCF);
					multiThreadRefCpgWriter = new SortingTcgaVCFWriter(refCpgWriter, MAXIMUM_CACHE_FOR_OUTPUT_VCF);
				}
			}
			
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        
    }

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.LocusWalker#map(org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker, org.broadinstitute.sting.gatk.contexts.ReferenceContext, org.broadinstitute.sting.gatk.contexts.AlignmentContext)
	 */
	@Override
	public VCFplusRef map(RefMetaDataTracker tracker, ReferenceContext ref,
			AlignmentContext context) {
		// TODO Auto-generated method stub
		
		if(tracker==null)
			return null;
		List<VariantContext> eval_bindings = tracker.getValues(eval);
		if(!eval_bindings.isEmpty()){
			
			
			VariantContext vc_eval = eval_bindings.get(0); 
			
			 return new VCFplusRef(ref, vc_eval);
		}
		
		return null;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduceInit()
	 */
	@Override
	public VCFCondition reduceInit() {
		// TODO Auto-generated method stub
		VCFCondition vcfCondition = new VCFCondition();
	//	initiateArray(vcfCondition.num);
	//	initiateArray(vcfCondition.methy);
		return vcfCondition;
	}

	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.Walker#reduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public VCFCondition reduce(VCFplusRef value, VCFCondition sum) {
		// TODO Auto-generated method stub
		
		if(value == null){
			return sum;
		}
		sum.nLociVisited++;
		if(islociFilterOut(value.vc)){
			sum.nLociFilterOut++;
			return sum;
		}
		//System.err.println(value.vc.getGenotype(0).getAttributes().keySet().toString());
		if(value.vc.hasAttribute(BisulfiteVCFConstants.CYTOSINE_TYPE) && !value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".").equalsIgnoreCase(".")){
			//System.err.println(value.vc.filtersWereApplied() + "\t" + value.vc.getGenotype(0).getGenotypeString() + "\t" + value.vc.getGenotype(0).getAttributes().keySet().toString());
			int numC=value.vc.getGenotype(0).getAttributeAsInt(BisulfiteVCFConstants.NUMBER_OF_C_KEY, -1);
			int numT=value.vc.getGenotype(0).getAttributeAsInt(BisulfiteVCFConstants.NUMBER_OF_T_KEY, -1);
			if(numC == -1 || numT == -1 || (numC==0 && numT==0))
				return sum;
			double methyValue = (double)numC/(double)(numC + numT);
			sum.nCytosineTotal++;
			sum.sumCytosineTotal += methyValue;
			//System.err.println(value.vc.getAttributeAsString(BisulfiteVCFConstants.CYTOSINE_TYPE, "."));
		//	if(value.vc.getAttributeAsString(BisulfiteVCFConstants.CYTOSINE_TYPE, ".").equalsIgnoreCase("[CG, CH, C]")){
			//	System.err.println(value.vc.getAttributeAsString(BisulfiteVCFConstants.CYTOSINE_TYPE, "."));
		//		return sum;
		//	}
				
			if(BisSNPUtils.isHomoC(value.vc)){
				if(BisSNPUtils.isRefCpg(value.refContext)){
					//String pattern = value.vc.getAttributeAsString(BisulfiteVCFConstants.CYTOSINE_TYPE, ".");
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					//System.err.println(value.vc.getAttributeAsString(BisulfiteVCFConstants.CYTOSINE_TYPE, "."));
					if(pattern.equalsIgnoreCase("CG")){
					//	sum.num[0][0][0]++;
					//	sum.methy[0][0][0] += methyValue;
						
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHom++;
						sum.sumCpgHom += methyValue;
						sum.nCpgHomRefCpg++;
						sum.sumCpgHomRefCpg += methyValue;
						
					}
					else if(pattern.equalsIgnoreCase("CH") ){
						
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHom++;
						sum.sumCphHom += methyValue;
						sum.nCphHomRefCpg++;
						sum.sumCphHomRefCpg += methyValue;
					}
					else if(pattern.equalsIgnoreCase(".") || pattern.equalsIgnoreCase("C")){
						
					}
					else{
						if(pattern.contains("C")){  //not only C but also contain some other base
						//if(BisSNPUtils.isHetCpg_at_G(value.vc)){
							sum.nCpgTotal++;
							sum.sumCpgTotal += methyValue;
							sum.nCpgHet++;
							sum.sumCpgHet += methyValue;
							sum.nCpgHetRefCpg++;
							sum.sumCpgHetRefCpg += methyValue;
							sum.nCpgHetAtG_RefCpg++;
							sum.sumCpgHetAtG_RefCpg += methyValue;
							if(outRefCphFile != null && outRefCpgFile != null)
								refCpgWriter.add(value.vc);
							//System.err.println(value.refContext.getLocus().toString() + "\t" + value.vc.getGenotype(0) + "\t" + numC + "\t" + numT + "\t" + sum.nCpgHetAtG_RefCpg + "\t" + sum.sumCpgHetAtG_RefCpg);
						}
					}
					
				}
				else if(BisSNPUtils.isRefCph(value.refContext)){
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if(pattern.equalsIgnoreCase("CG")){
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHom++;
						sum.sumCpgHom += methyValue;
						sum.nCpgHomRefCph++;
						sum.sumCpgHomRefCph += methyValue;
					//	System.err.println(value.refContext.getLocus().toString() + "\t" + value.vc.getGenotype(0) + "\t" + numC + "\t" + numT + "\t" + sum.nCpgHomRefCph + "\t" + sum.sumCpgHomRefCph);
					}
					else if(pattern.equalsIgnoreCase("CH") ){
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHom++;
						sum.sumCphHom += methyValue;
						sum.nCphHomRefCph++;
						sum.sumCphHomRefCph += methyValue;
					}
					else if(pattern.equalsIgnoreCase(".") || pattern.equalsIgnoreCase("C")){
						
					}
					else{
						if(pattern.contains("C")){
							sum.nCpgTotal++;
							sum.sumCpgTotal += methyValue;
							sum.nCpgHet++;
							sum.sumCpgHet += methyValue;
							sum.nCpgHetRefCph++;
							sum.sumCpgHetRefCph += methyValue;
							sum.nCpgHetAtG_RefCph++;
							sum.sumCpgHetAtG_RefCph += methyValue;
							if(outRefCphFile != null && outRefCpgFile != null)
								refCphWriter.add(value.vc);
						}
					}
				}
				else{
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if(pattern.equalsIgnoreCase("CG")){
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHom++;
						sum.sumCpgHom += methyValue;
						sum.nCpgHomRefNotC++;
						sum.sumCpgHomRefNotC += methyValue;
					}
					else if(pattern.equalsIgnoreCase("CH") ){
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHom++;
						sum.sumCphHom += methyValue;
						sum.nCphHomRefNotC++;
						sum.sumCphHomRefNotC += methyValue;
					}
					else if(pattern.equalsIgnoreCase(".") || pattern.equalsIgnoreCase("C")){
						
					}
					else{
						if(pattern.contains("C")){
							sum.nCpgTotal++;
							sum.sumCpgTotal += methyValue;
							sum.nCpgHet++;
							sum.sumCpgHet += methyValue;
							sum.nCpgHetRefNotC++;
							sum.sumCpgHetRefNotC += methyValue;
							sum.nCpgHetAtG_RefNotC++;
							sum.sumCpgHetAtG_RefNotC += methyValue;
						}
					}
				}
			}
			else{//exclude C/T SNP at cytosine sites!!
				if(BisSNPUtils.isRefCpg(value.refContext)){
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if(pattern.contains("G") && !pattern.contains("Y")){
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHet++;
						sum.sumCpgHet += methyValue;
						sum.nCpgHetRefCpg++;
						sum.sumCpgHetRefCpg += methyValue;
						sum.nCpgHetAtC_RefCpg++;
						sum.sumCpgHetAtC_RefCpg += methyValue;
					}
					else if(pattern.contains("H") && !pattern.contains("Y")){
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHet++;
						sum.sumCphHet += methyValue;
						sum.nCphHetRefCpg++;
						sum.sumCphHetRefCpg += methyValue;
						sum.nCphHetAtC_RefCpg++;
						sum.sumCphHetAtC_RefCpg += methyValue;
					}
					
				}
				else if(BisSNPUtils.isRefCph(value.refContext)){
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if(pattern.contains("G") && !pattern.contains("Y")){
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHet++;
						sum.sumCpgHet += methyValue;
						sum.nCpgHetRefCph++;
						sum.sumCpgHetRefCph += methyValue;
						sum.nCpgHetAtC_RefCph++;
						sum.sumCpgHetAtC_RefCph += methyValue;
					}
					else if(pattern.contains("H") && !pattern.contains("Y")){
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHet++;
						sum.sumCphHet += methyValue;
						sum.nCphHetRefCph++;
						sum.sumCphHetRefCph += methyValue;
						sum.nCphHetAtC_RefCph++;
						sum.sumCphHetAtC_RefCph += methyValue;
					}
				}
				else{
					String pattern = value.vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.BEST_C_PATTERN, ".");
					if(pattern.contains("G") && !pattern.contains("Y")){
						sum.nCpgTotal++;
						sum.sumCpgTotal += methyValue;
						sum.nCpgHet++;
						sum.sumCpgHet += methyValue;
						sum.nCpgHetRefNotC++;
						sum.sumCpgHetRefNotC += methyValue;
						sum.nCpgHetAtC_RefNotC++;
						sum.sumCpgHetAtC_RefNotC += methyValue;
					}
					else if(pattern.contains("H") && !pattern.contains("Y")){
						sum.nCphTotal++;
						sum.sumCphTotal += methyValue;
						sum.nCphHet++;
						sum.sumCphHet += methyValue;
						sum.nCphHetRefNotC++;
						sum.sumCphHetRefNotC += methyValue;
						sum.nCphHetAtC_RefNotC++;
						sum.sumCphHetAtC_RefNotC += methyValue;
					}
				}
			}
		}
		
		return sum;
	}
	
	
	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 */
	@Override
	public VCFCondition treeReduce(VCFCondition lhs, VCFCondition rhs) {
		// TODO Auto-generated method stub
		lhs.nLociVisited += rhs.nLociVisited;
		lhs.nLociFilterOut += rhs.nLociFilterOut;
        lhs.nCytosineTotal += rhs.nCytosineTotal;
        
		lhs.nCpgTotal += rhs.nCpgTotal;
		lhs.nCpgHom += rhs.nCpgHom;
		lhs.nCpgHomRefCpg += rhs.nCpgHomRefCpg;
		lhs.nCpgHomRefCph += rhs.nCpgHomRefCph;
		lhs.nCpgHomRefNotC += rhs.nCpgHomRefNotC;
		lhs.nCpgHet += rhs.nCpgHet;
		lhs.nCpgHetRefCpg += rhs.nCpgHetRefCpg;
		lhs.nCpgHetRefCph += rhs.nCpgHetRefCph;
		
		lhs.nCpgHetAtC_RefCpg += rhs.nCpgHetAtC_RefCpg;
		lhs.nCpgHetAtG_RefCpg += rhs.nCpgHetAtG_RefCpg;
		lhs.nCpgHetAtC_RefCph += rhs.nCpgHetAtC_RefCph;
		lhs.nCpgHetAtG_RefCph += rhs.nCpgHetAtG_RefCph;
		lhs.nCpgHetAtC_RefNotC += rhs.nCpgHetAtC_RefNotC;
		lhs.nCpgHetAtG_RefNotC += rhs.nCpgHetAtG_RefNotC;
		
		lhs.nCphTotal += rhs.nCphTotal;
		lhs.nCphHom += rhs.nCphHom;
		lhs.nCphHomRefCpg += rhs.nCphHomRefCpg;
		lhs.nCphHomRefCph += rhs.nCphHomRefCph;
		lhs.nCphHomRefNotC += rhs.nCphHomRefNotC;
		lhs.nCphHet += rhs.nCphHet;
		lhs.nCphHetRefCpg += rhs.nCphHetRefCpg;
		lhs.nCphHetRefCph += rhs.nCphHetRefCph;
		
		lhs.nCphHetAtC_RefCpg += rhs.nCphHetAtC_RefCpg;
		lhs.nCphHetAtC_RefCph += rhs.nCphHetAtC_RefCph;
		lhs.nCphHetAtC_RefNotC += rhs.nCphHetAtC_RefNotC;
		lhs.nCphHetAtC_RefNotC += rhs.nCphHetAtC_RefNotC;
		
		lhs.sumCytosineTotal += rhs.sumCytosineTotal;
        
		lhs.sumCpgTotal += rhs.sumCpgTotal;
		lhs.sumCpgHom += rhs.sumCpgHom;
		lhs.sumCpgHomRefCpg += rhs.sumCpgHomRefCpg;
		lhs.sumCpgHomRefCph += rhs.sumCpgHomRefCph;
		lhs.sumCpgHomRefNotC += rhs.sumCpgHomRefNotC;
		lhs.sumCpgHet += rhs.sumCpgHet;
		lhs.sumCpgHetRefCpg += rhs.sumCpgHetRefCpg;
		lhs.sumCpgHetRefCph += rhs.sumCpgHetRefCph;
		
		lhs.sumCpgHetAtC_RefCpg += rhs.sumCpgHetAtC_RefCpg;
		lhs.sumCpgHetAtG_RefCpg += rhs.sumCpgHetAtG_RefCpg;
		lhs.sumCpgHetAtC_RefCph += rhs.sumCpgHetAtC_RefCph;
		lhs.sumCpgHetAtG_RefCph += rhs.sumCpgHetAtG_RefCph;
		lhs.sumCpgHetAtC_RefNotC += rhs.sumCpgHetAtC_RefNotC;
		lhs.sumCpgHetAtG_RefNotC += rhs.sumCpgHetAtG_RefNotC;
		
		lhs.sumCphTotal += rhs.sumCphTotal;
		lhs.sumCphHom += rhs.sumCphHom;
		lhs.sumCphHomRefCpg += rhs.sumCphHomRefCpg;
		lhs.sumCphHomRefCph += rhs.sumCphHomRefCph;
		lhs.sumCphHomRefNotC += rhs.sumCphHomRefNotC;
		lhs.sumCphHet += rhs.sumCphHet;
		lhs.sumCphHetRefCpg += rhs.sumCphHetRefCpg;
		lhs.sumCphHetRefCph += rhs.sumCphHetRefCph;
		
		lhs.sumCphHetAtC_RefCpg += rhs.sumCphHetAtC_RefCpg;
		lhs.sumCphHetAtC_RefCph += rhs.sumCphHetAtC_RefCph;
		lhs.sumCphHetAtC_RefNotC += rhs.sumCphHetAtC_RefNotC;
		lhs.sumCphHetAtC_RefNotC += rhs.sumCphHetAtC_RefNotC;
		
		return lhs;
	}

	public void onTraversalDone(VCFCondition result) {
		//num:     c, cpg, cph, cpgHom, cpgHomRefCpg, cpgHomRef
		//methyAve:
		//all
		//writer.println("nLociVisited\t" + result.nLociVisited);
		writer.println("nCytosineTotal\t" + result.nCytosineTotal + "\tsumCytosineTotal\t" + 100*result.sumCytosineTotal/result.nCytosineTotal);
		writer.println("nCpgTotal\t" + result.nCpgTotal + "\tsumCpgTotal\t" + 100*result.sumCpgTotal/result.nCpgTotal);
		writer.println("nCphTotal\t" + result.nCphTotal + "\tsumCphTotal\t" + 100*result.sumCphTotal/result.nCphTotal);
		
		//hom cytosine:
		writer.println("nCpgHom\t" + result.nCpgHom + "\t" + "sumCpgHom\t" + 100*result.sumCpgHom/result.nCpgHom);
		writer.println("nCpgHomRefCpg\t" + result.nCpgHomRefCpg + "\t" + "sumCpgHomRefCpg\t" + 100*result.sumCpgHomRefCpg/result.nCpgHomRefCpg);
		writer.println("nCpgHomRefCph\t" + result.nCpgHomRefCph + "\t" + "sumCpgHomRefCph\t" + 100*result.sumCpgHomRefCph/result.nCpgHomRefCph);
		writer.println("nCpgHomRefNotC\t" + result.nCpgHomRefNotC + "\t" + "sumCpgHomRefNotC\t" + 100*result.sumCpgHomRefNotC/result.nCpgHomRefNotC);
		
		writer.println("nCphHom\t" + result.nCphHom + "\t" + "sumCphHom\t" + 100*result.sumCphHom/result.nCphHom);
		writer.println("nCphHomRefCpg\t" + result.nCphHomRefCpg + "\t" + "sumCphHomRefCpg\t" + 100*result.sumCphHomRefCpg/result.nCphHomRefCpg);
		writer.println("nCphHomRefCph\t" + result.nCphHomRefCph + "\t" + "sumCphHomRefCph\t" + 100*result.sumCphHomRefCph/result.nCphHomRefCph);
		writer.println("nCphHomRefNotC\t" + result.nCphHomRefNotC + "\t" + "sumCphHomRefNotC\t" + 100*result.sumCphHomRefNotC/result.nCphHomRefNotC);
		
		//het cytosines:
		writer.println("nCpgHet\t" + result.nCpgHet + "\t" + "sumCpgHet\t" + 100*result.sumCpgHet/result.nCpgHet);
		writer.println("nCpgHetRefCpg\t" + result.nCpgHetRefCpg + "\t" + "sumCpgHetRefCpg\t" + 100*result.sumCpgHetRefCpg/result.nCpgHetRefCpg);
		writer.println("nCpgHetRefCph\t" + result.nCpgHetRefCph + "\t" + "sumCpgHetRefCph\t" + 100*result.sumCpgHetRefCph/result.nCpgHetRefCph);
		writer.println("nCpgHetAtC_RefCpg\t" + result.nCpgHetAtC_RefCpg + "\t" + "sumCpgHetAtC_RefCpg\t" + 100*result.sumCpgHetAtC_RefCpg/result.nCpgHetAtC_RefCpg);
		writer.println("nCpgHetAtC_RefCph\t" + result.nCpgHetAtC_RefCph + "\t" + "sumCpgHetAtC_RefCph\t" + 100*result.sumCpgHetAtC_RefCph/result.nCpgHetAtC_RefCph);
		writer.println("nCpgHetAtC_RefNotC\t" + result.nCpgHetAtC_RefNotC + "\t" + "sumCpgHetAtC_RefNotC\t" + 100*result.sumCpgHetAtC_RefNotC/result.nCpgHetAtC_RefNotC);
		writer.println("nCpgHetAtG_RefCpg\t" + result.nCpgHetAtG_RefCpg + "\t" + "sumCpgHetAtG_RefCpg\t" + 100*result.sumCpgHetAtG_RefCpg/result.nCpgHetAtG_RefCpg);
		writer.println("nCpgHetAtG_RefCph\t" + result.nCpgHetAtG_RefCph + "\t" + "sumCpgHetAtG_RefCph\t" + 100*result.sumCpgHetAtG_RefCph/result.nCpgHetAtG_RefCph);
		writer.println("nCpgHetAtG_RefNotC\t" + result.nCpgHetAtG_RefNotC + "\t" + "sumCpgHetAtG_RefNotC\t" + 100*result.sumCpgHetAtG_RefNotC/result.nCpgHetAtG_RefNotC);
		
		writer.println("nCphHet\t" + result.nCphHet + "\t" + "sumCphHet\t" + 100*result.sumCphHet/result.nCphHet);
		writer.println("nCphHetRefCpg\t" + result.nCphHetRefCpg + "\t" + "sumCphHetRefCpg\t" + 100*result.sumCphHetRefCpg/result.nCphHetRefCpg);
		writer.println("nCphHetRefCph\t" + result.nCphHetRefCph + "\t" + "sumCphHetRefCph\t" + 100*result.sumCphHetRefCph/result.nCphHetRefCph);
		writer.println("nCphHetAtC_RefCpg\t" + result.nCphHetAtC_RefCpg + "\t" + "sumCphHetAtC_RefCpg\t" + 100*result.sumCphHetAtC_RefCpg/result.nCphHetAtC_RefCpg);
		writer.println("nCphHetAtC_RefCph\t" + result.nCphHetAtC_RefCph + "\t" + "sumCphHetAtC_RefCph\t" + 100*result.sumCphHetAtC_RefCph/result.nCphHetAtC_RefCph);
		writer.println("nCphHetAtC_RefNotC\t" + result.nCphHetAtC_RefNotC + "\t" + "sumCphHetAtC_RefNotC\t" + 100*result.sumCphHetAtC_RefNotC/result.nCphHetAtC_RefNotC);
		
		writer.close();
		if(outRefCphFile != null && outRefCpgFile != null){
			refCphWriter.close();
			refCpgWriter.close();
			if(getToolkit().getArguments().numberOfThreads>1){
				multiThreadRefCphWriter.close();
				multiThreadRefCpgWriter.close();
			}
		}
		
		logger.info("Finished!");
		logger.info(String.format("Visited Loci: %d,\t Filtered Loci: %d", result.nLociVisited, result.nLociFilterOut));
		
	}

	public class VCFplusRef{
		public ReferenceContext refContext;
		public VariantContext vc;
		public VCFplusRef(ReferenceContext ref, VariantContext vc_eval){
			refContext = ref;
			vc = vc_eval;
		}
		
	}
	
	private void initiateArray(long[][][] arrays){
		for(int i=0; i<arrays.length; i++){
			for(int ii=0; i<arrays[0].length; ii++){
				for(int iii=0; i<arrays[0][0].length; iii++){
					arrays[i][ii][iii]=0;
				}
			}
		}
	}
	
	private void initiateArray(double[][][] arrays){
		for(int i=0; i<arrays.length; i++){
			for(int ii=0; i<arrays[0].length; ii++){
				for(int iii=0; i<arrays[0][0].length; iii++){
					arrays[i][ii][iii]=0.0;
				}
			}
		}
	}
	
	private boolean islociFilterOut(VariantContext vc){
		if(vc.getPhredScaledQual() <= qual)
			return true;
		if(vc.hasAttribute(VCFConstants.DEPTH_KEY) && vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > maxCov){
			return true;
		}
		if(vc.hasAttribute(VCFConstants.STRAND_BIAS_KEY) && vc.getAttributeAsDouble(VCFConstants.STRAND_BIAS_KEY, 0.0) > sb){
			return true;
		}
		if(vc.hasAttribute(VCFConstants.MAPPING_QUALITY_ZERO_KEY)){
			if(vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > 40 && vc.getAttributeAsDouble(VCFConstants.MAPPING_QUALITY_ZERO_KEY, 0)/vc.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0) > mq0)
				return true;
		}
		if(vc.hasAttribute(BisulfiteVCFConstants.QUAL_BY_DEPTH) && vc.getAttributeAsDouble(BisulfiteVCFConstants.QUAL_BY_DEPTH, 10000.0) < qd){
			return true;
		}
		if(vc.getGenotype(0).hasAttribute(BisulfiteVCFConstants.C_STATUS)){
			String cytosineStatus = vc.getGenotype(0).getAttributeAsString(BisulfiteVCFConstants.C_STATUS, ".");
			if(!cytosineStatus.equalsIgnoreCase(".")){
				String[] tmp = cytosineStatus.split(",");
				int numC = Integer.parseInt(tmp[0]);
				int numT = Integer.parseInt(tmp[1]);
				//System.err.println(cytosineStatus + "\t" + tmp + "\t" + numC+ "\t" + numT);
				if(numC + numT < minCT){
					return true;
				}
			}
			
		}
		return false;
	}
}
