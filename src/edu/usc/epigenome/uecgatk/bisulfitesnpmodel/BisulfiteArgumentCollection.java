package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;

public class BisulfiteArgumentCollection extends UnifiedArgumentCollection {
	@Argument(fullName = "sequencing_mode", shortName = "sm", doc = "Bisulfite mode: BM, GNOMe-seq mode: GM, Normal sequencing mode: NM", required = false)
    public NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel sequencingMode = NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel.BM;
	
	@Argument(fullName = "auto_estimate_cpg_methylation", shortName = "aecpg", doc = "the first run would be to run auto_estimate_cpg methylation status", required = false)
    public boolean autoEstimateCpg = true;
	
	@Argument(fullName = "auto_estimate_chg_methylation", shortName = "aechg", doc = "the first run would be to run auto_estimate_chg methylation status", required = false)
    public boolean autoEstimateChg = true;
	
	@Argument(fullName = "auto_estimate_chh_methylation", shortName = "aechh", doc = "the first run would be to run auto_estimate_chh methylation status", required = false)
    public boolean autoEstimateChh = true;
	
	@Argument(fullName = "auto_estimate_gch_methylation", shortName = "aegch", doc = "the first run would be to run auto_estimate_gch methylation status", required = false)
    public boolean autoEstimateGch = true;
	
	@Argument(fullName = "auto_estimate_gcg_methylation", shortName = "aegcg", doc = "the first run would be to run auto_estimate_gcg methylation status", required = false)
    public boolean autoEstimateGcg = true;
	
	@Argument(fullName = "auto_estimate_hcg_methylation", shortName = "aehcg", doc = "the first run would be to run auto_estimate_hcg methylation status", required = false)
    public boolean autoEstimateHcg = true;
	
	@Argument(fullName = "auto_estimate_other_cytosine_methylation", shortName = "aeoc", doc = "the first run would be to run auto_estimate_other_cytosine_methylation status, you need to provide cytosine type by such format(GCAA is ctosine type, 2 means cytosine is in 2nd base): -aoec GCAA:2;GGGCA:4", required = false)
    //public String autoEstimateOtherCytosine = "GCAA-2;GGGCA-4";
	public String autoEstimateOtherCytosine = "";
	
	@Argument(fullName = "force_cpg_methylation", shortName = "fcpg", doc = "force the cpg methylation status", required = false)
    public double forceCpg = 0.50;
	
	@Argument(fullName = "force_chg_methylation", shortName = "fchg", doc = "force the chg methylation status", required = false)
    public double forceChg = 0.50;
	
	@Argument(fullName = "force_chh_methylation", shortName = "fchh", doc = "force the chh methylation status", required = false)
    public double forceChh = 0.50;
	
	@Argument(fullName = "force_gch_methylation", shortName = "fgch", doc = "force the gch methylation status", required = false)
    public double forceGch = 0.50;
	
	@Argument(fullName = "force_gcg_methylation", shortName = "fgcg", doc = "force the gcg methylation status", required = false)
    public double forceGcg = 0.50;
	
	@Argument(fullName = "force_hcg_methylation", shortName = "fhcg", doc = "force the hcg methylation status", required = false)
    public double forceHcg = 0.50;
	
	@Argument(fullName = "force_other_cytosine_methylation", shortName = "foc", doc = "force the other_cytosine_methylation status, you need to provide cytosine type by such format(GCAA is ctosine type, 2 means cytosine is in 2nd base, 0.75 means methylation level): -aoec GCAA:2:0.75;GGGCA:3:0.33", required = false)
    //public String forceOtherCytosine = "GCAA-2:0.75;GGGCA-4:0.33";
	public String forceOtherCytosine = "";
	
//need to improve..
	@Argument(fullName = "log_likelihood_ratio_for_cytosine_type", shortName = "cTypeThreshold", doc = "phred scale likelihood ratio of threshold to be this cytosine type but not other cytosine, default is 20, means 100 times more than the other type of cytosine", required = false)
    public double cTypeThreshold = 10;
	
	//@Argument(fullName = "Cytosine_Type", shortName = "ct", doc = "Cytosine type, CG, CHH, CHG or GCH....for test only (format should be -ct CG-0:0.75;CHH-0:0.01... add the cytosine type, cytosine position in your string and their genome wide methylation value you estimate )", required = false)
    //public String cytosineType = "CGA-0:0.7314;GCA-1:0.01";
	//public String cytosineType = null;
	
	@Argument(fullName = "test_location", shortName = "loc", doc = ".for test only", required = false)
    public long testLocus = -1;
	
	@Argument(fullName = "bisulfite_conversion_rate", shortName = "bsRate", doc = "bisulfite_conversion_rate .for test only", required = false)
    public double bsRate = 0.9975;
	
	@Argument(fullName = "over_conversion_rate", shortName = "overRate", doc = "cytosine_over_conversion_rate .for test only", required = false)
    public double overRate = 0;
	
	//@Argument(fullName = "CpG_Methylation_rate_in_CGI", shortName = "CpgMethyCGI", doc = "CpG_Methylation_rate_in_CGI .for test only", required = false)
    //public double CpgMethyCGI = 0;
	
//	@Argument(fullName = "CpG_Methylation_rate_not_in_CGI", shortName = "CpgMethyNonCGI", doc = "CpG_Methylation_rate_not_in_CGI .for test only", required = false)
 //   public double CpgMethyNonCGI = 0;
	
//	@Argument(fullName = "CpH_Methylation_rate", shortName = "CphMethy", doc = "CpH_Methylation_rate .for test only", required = false)
 //   public double CphMethy = 0;
	
	@Argument(fullName = "validateDbsnphet", shortName = "vdh", doc = "validateDbsnphet .for test only", required = false)
    public double validateDbsnpHet = 0.1;
	
	@Argument(fullName = "novelDbsnpHet", shortName = "ndh", doc = "novelDbsnpHet .for test only", required = false)
    public double novelDbsnpHet = 0.02;
	
	
	public BisulfiteArgumentCollection clone() {
		BisulfiteArgumentCollection bac = new BisulfiteArgumentCollection();
		bac.GLmodel = GLmodel;
        bac.heterozygosity = heterozygosity;
        bac.PCR_error = PCR_error;
        bac.GenotypingMode = GenotypingMode;
        bac.OutputMode = OutputMode;
        bac.NO_SLOD = NO_SLOD;
        bac.ASSUME_SINGLE_SAMPLE = ASSUME_SINGLE_SAMPLE;
        bac.STANDARD_CONFIDENCE_FOR_CALLING = STANDARD_CONFIDENCE_FOR_CALLING;
        bac.STANDARD_CONFIDENCE_FOR_EMITTING = STANDARD_CONFIDENCE_FOR_EMITTING;
        bac.MIN_BASE_QUALTY_SCORE = MIN_BASE_QUALTY_SCORE;
        bac.MIN_MAPPING_QUALTY_SCORE = MIN_MAPPING_QUALTY_SCORE;
        bac.MAX_MISMATCHES = MAX_MISMATCHES;
        bac.USE_BADLY_MATED_READS = USE_BADLY_MATED_READS;
        bac.MAX_DELETION_FRACTION = MAX_DELETION_FRACTION;
        bac.MIN_INDEL_COUNT_FOR_GENOTYPING = MIN_INDEL_COUNT_FOR_GENOTYPING;
        bac.INDEL_HETEROZYGOSITY = INDEL_HETEROZYGOSITY;
        bac.INSERTION_START_PROBABILITY = INSERTION_START_PROBABILITY;
        bac.INSERTION_END_PROBABILITY = INSERTION_END_PROBABILITY;
        bac.ALPHA_DELETION_PROBABILITY = ALPHA_DELETION_PROBABILITY;
        bac.sequencingMode = sequencingMode;
        bac.autoEstimateChg = autoEstimateChg;
        bac.autoEstimateChh = autoEstimateChh;
        bac.autoEstimateCpg = autoEstimateCpg;
        bac.autoEstimateGch = autoEstimateGch;
        bac.autoEstimateGcg = autoEstimateGcg;
        bac.autoEstimateHcg = autoEstimateHcg;
        bac.forceChg = forceChg;
        bac.forceCpg = forceCpg;
        bac.forceChh = forceChh;
        bac.forceGch = forceGch;
        bac.forceGcg = forceGcg;
        bac.forceHcg = forceHcg;
        
       // bac.cytosineType = cytosineType;
        bac.testLocus = testLocus;
        bac.bsRate = bsRate;
        bac.overRate = overRate;
        //bac.CpgMethyCGI = CpgMethyCGI;
       // bac.CpgMethyNonCGI = CpgMethyNonCGI;
       // bac.CphMethy = CphMethy;
        bac.validateDbsnpHet = validateDbsnpHet;
        bac.novelDbsnpHet = novelDbsnpHet;
        
        return bac;
    }

}
