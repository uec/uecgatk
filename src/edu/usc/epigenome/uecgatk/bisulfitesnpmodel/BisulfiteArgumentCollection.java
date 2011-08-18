package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;

public class BisulfiteArgumentCollection extends UnifiedArgumentCollection {
	@Argument(fullName = "sequencing_mode", shortName = "sm", doc = "Bisulfite mode: BM, GNOMe-seq mode: GM, Normal sequencing mode: NM", required = false)
    public NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel sequencingMode = NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel.BM;
	
	@Argument(fullName = "paired_end_mode", shortName = "pem", doc = "work in paired end mode", required = false)
    public boolean pairedEndMode = false;
	
	@Argument(fullName = "auto_estimate_cpg_methylation", shortName = "aecpg", doc = "the first run would be to run auto_estimate_cpg methylation status", required = false)
    public boolean autoEstimateCpg = true;
	
	@Argument(fullName = "auto_estimate_cph_methylation", shortName = "aecph", doc = "the first run would be to run auto_estimate_cph methylation status", required = false)
    public boolean autoEstimateCph = true;
	
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
	
	@Argument(fullName = "auto_estimate_other_cytosine_methylation", shortName = "aeoc", doc = "the first run would be to run auto_estimate_other_cytosine_methylation status, you need to provide cytosine type by such format: -aoec GCAA-2:0.5;GGGCA-4:0.5 ((GCAA is ctosine type, 2 means cytosine is in 2nd base, 0.5 means intial methylation status))", required = false)
    //example: "GCAA-2:0.5;GGGCA-4:0.5";
	public String autoEstimateOtherCytosine = "";
	
	@Argument(fullName = "force_cpg_methylation", shortName = "fcpg", doc = "force the cpg methylation status", required = false)
    public double forceCpg = 0.50;
	
	@Argument(fullName = "force_cph_methylation", shortName = "fcph", doc = "force the cph methylation status", required = false)
    public double forceCph = 0.50;
	
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
	
	@Argument(fullName = "force_other_cytosine_methylation", shortName = "foc", doc = "force the other_cytosine_methylation status, you need to provide cytosine type by such format: -aoec GCAA-2:0.75;GGGCA-3:0.33 (GCAA is ctosine type, 2 means cytosine is in 2nd base, 0.75 means methylation level)", required = false)
    //example String forceOtherCytosine = "GCAA-2:0.75;GGGCA-4:0.33";
	public String forceOtherCytosine = "";
	

	@Argument(fullName = "log_likelihood_ratio_for_cytosine_type", shortName = "cTypeThreshold", doc = "phred scale likelihood ratio of to be this cytosine pattern but not other cytosines in the first iteration for two-iteration mode (the real criteria is cTypeThreshold + stand_call_conf), default is 20, if stand_call_conf is 0, means 10^((20+0)/10) = 100 times more likihood than the other type of cytosine, only used in the first iteration", required = false)
    public double cTypeThreshold = 20;
	
	@Argument(fullName = "test_location", shortName = "loc", doc = "for debug only, output the detail information in the location", required = false)
    public long testLocus = -1;
	
	@Argument(fullName = "bisulfite_conversion_rate", shortName = "bsRate", doc = "bisulfite conversion rate", required = false)
    public double bsRate = 0.9975;
	
	@Argument(fullName = "over_conversion_rate", shortName = "overRate", doc = "cytosine over conversion rate. it is often 0", required = false)
    public double overRate = 0;
	
	@Argument(fullName = "validateDbsnphet", shortName = "vdh", doc = "heterozygous SNP rate when the loci is discovered as SNP in dbSNP and is validated, the default value is human genome", required = false)
    public double validateDbsnpHet = 0.1;
	
	@Argument(fullName = "novelDbsnpHet", shortName = "ndh", doc = "heterozygous SNP rate when the loci is discovered as SNP in dbSNP and but not validated, the default value is human genome", required = false)
    public double novelDbsnpHet = 0.02;
	
	@Argument(fullName = "reference_genome_error", shortName = "rge", doc = "reference genome error, the default value is human genome, in hg17, it is less than 1e-4; in SOAPsnp, it is 1e-5; in GATK it is 1e-6. so is it because hg18 error rate is 1e-5, and hg19 is 1e-6? can't find any reference about it..", required = false)
    public double referenceGenomeErr = 1e-6;
	
	@Argument(fullName = "tcga_format_vcf", shortName = "tcga", doc = "output TCGA specific VCF format or not, not used yet, in test", required = false)
    public boolean tcga = false;
	
    @Argument(fullName = "output_genotype", shortName = "out_genotype", doc = "Should we output confident genotypes (i.e. including ref calls),just the variants, just homozygous CpG or just homozygous Cytosines?", required = false)
    public BisulfiteGenotyperEngine.OUTPUT_MODE OutputMode = BisulfiteGenotyperEngine.OUTPUT_MODE.EMIT_ALL_CONFIDENT_SITES;
	
	
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
        bac.pairedEndMode = pairedEndMode;
        bac.autoEstimateChg = autoEstimateChg;
        bac.autoEstimateChh = autoEstimateChh;
        bac.autoEstimateCpg = autoEstimateCpg;
        bac.autoEstimateCph = autoEstimateCph;
        bac.autoEstimateGch = autoEstimateGch;
        bac.autoEstimateGcg = autoEstimateGcg;
        bac.autoEstimateHcg = autoEstimateHcg;
        bac.forceChg = forceChg;
        bac.forceCpg = forceCpg;
        bac.forceCph = forceCph;
        bac.forceChh = forceChh;
        bac.forceGch = forceGch;
        bac.forceGcg = forceGcg;
        bac.forceHcg = forceHcg;
        bac.autoEstimateOtherCytosine = autoEstimateOtherCytosine;
        bac.forceOtherCytosine = forceOtherCytosine;
        
        bac.cTypeThreshold = cTypeThreshold;
        bac.testLocus = testLocus;
        bac.bsRate = bsRate;
        bac.overRate = overRate;
        bac.validateDbsnpHet = validateDbsnpHet;
        bac.novelDbsnpHet = novelDbsnpHet;
        bac.referenceGenomeErr = referenceGenomeErr;
        bac.heterozygosity = heterozygosity;
    
        return bac;
    }

}
