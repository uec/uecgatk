package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.lang.reflect.Field;

import org.broadinstitute.sting.commandline.Argument;
import edu.usc.epigenome.uecgatk.bisulfitesnpmodel.NonRefDependSNPGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;

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

public class BisulfiteArgumentCollection extends UnifiedArgumentCollection {
	@Argument(fullName = "sequencing_mode", shortName = "sm", doc = "Bisulfite mode: BM, GNOMe-seq mode: GM, Normal sequencing mode: NM", required = false)
    public NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel sequencingMode = NonRefDependSNPGenotypeLikelihoodsCalculationModel.MethylSNPModel.BM;
	
	@Argument(fullName = "paired_end_mode", shortName = "pem", doc = "work in paired end mode", required = false)
    public boolean pairedEndMode = false;
	
	@Argument(fullName = "bisulfite_conversion_only_on_one_strand", shortName = "bcm", doc = "true: Illumina protocol which is often used, only bisulfite conversion strand is kept (Lister protocol, sequence 2 forward strands only); false: Cokus protocol, sequence all 4 bisulfite converted strands", required = false)
    public boolean bisulfiteConversionModeOnestrand = true;
	
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
	
	@Argument(fullName = "auto_estimate_hch_methylation", shortName = "aehch", doc = "the first run would be to run auto_estimate_hch methylation status", required = false)
    public boolean autoEstimateHch = true;
	
	@Argument(fullName = "auto_estimate_wch_methylation", shortName = "aewch", doc = "the first run would be to run auto_estimate_wch methylation status", required = false)
    public boolean autoEstimateWch = true;
	
	@Argument(fullName = "auto_estimate_gcg_methylation", shortName = "aegcg", doc = "the first run would be to run auto_estimate_gcg methylation status", required = false)
    public boolean autoEstimateGcg = true;
	
	@Argument(fullName = "auto_estimate_hcg_methylation", shortName = "aehcg", doc = "the first run would be to run auto_estimate_hcg methylation status", required = false)
    public boolean autoEstimateHcg = true;
	
	@Argument(fullName = "auto_estimate_wcg_methylation", shortName = "aewcg", doc = "the first run would be to run auto_estimate_wcg methylation status", required = false)
    public boolean autoEstimateWcg = true;
	
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
	
	@Argument(fullName = "force_hch_methylation", shortName = "fhch", doc = "force the hch methylation status", required = false)
    public double forceHch = 0.50;
	
	@Argument(fullName = "force_wch_methylation", shortName = "fwch", doc = "force the wch methylation status", required = false)
    public double forceWch = 0.50;
	
	@Argument(fullName = "force_gcg_methylation", shortName = "fgcg", doc = "force the gcg methylation status", required = false)
    public double forceGcg = 0.50;
	
	@Argument(fullName = "force_hcg_methylation", shortName = "fhcg", doc = "force the hcg methylation status", required = false)
    public double forceHcg = 0.50;
	
	@Argument(fullName = "force_wcg_methylation", shortName = "fwcg", doc = "force the wcg methylation status", required = false)
    public double forceWcg = 0.50;
	
	@Argument(fullName = "force_other_cytosine_methylation", shortName = "foc", doc = "force the other_cytosine_methylation status, you need to provide cytosine type by such format: -aoec GCAA-2:0.75;GGGCA-3:0.33 (GCAA is ctosine type, 2 means cytosine is in 2nd base, 0.75 means methylation level)", required = false)
    //example String forceOtherCytosine = "GCAA-2:0.75;GGGCA-4:0.33";
	public String forceOtherCytosine = "";
	

	@Argument(fullName = "log_likelihood_ratio_for_cytosine_type", shortName = "cTypeThreshold", doc = "phred scale likelihood ratio of to be this cytosine pattern but not other cytosines in the first iteration for two-iteration mode (the real criteria is cTypeThreshold + stand_call_conf), default is 20, if stand_call_conf is 0, means 10^((20+0)/10) = 100 times more likihood than the other type of cytosine, only used in the first iteration", required = false)
    public double cTypeThreshold = 20;
	
	@Argument(fullName = "test_location", shortName = "loc", doc = "for debug only, output the detail information in the location", required = false)
    public long testLocus = -1;
	
	@Argument(fullName = "minmum_cytosine_converted", shortName = "minConv", doc = "disregard first few cytosines in the reads which may come from uncomplete bisulfite conversion in the first few cytosines of the reads", required = false)
    public short minConv = 0;
	
	@Argument(fullName = "bisulfite_conversion_rate", shortName = "bsRate", doc = "bisulfite conversion rate", required = false)
    public double bsRate = 0.9975;
	
	@Argument(fullName = "over_conversion_rate", shortName = "overRate", doc = "cytosine over conversion rate. it is often 0", required = false)
    public double overRate = 0;
	
	@Argument(fullName = "validateDbsnphet", shortName = "vdh", doc = "heterozygous SNP rate when the loci is discovered as SNP in dbSNP and is validated, the default value is human genome", required = false)
    public double validateDbsnpHet = 0.1;
	
	@Argument(fullName = "novelDbsnpHet", shortName = "ndh", doc = "heterozygous SNP rate when the loci is discovered as SNP in dbSNP and but not validated, the default value is human genome", required = false)
    public double novelDbsnpHet = 0.02;
	
	@Argument(fullName = "reference_genome_error", shortName = "rge", doc = "Reference genome error, the default value is human genome, in hg16 it is 99.99% accurate,  in hg17/hg18/hg19, it is less than 1e-4 (USCS genome browser described); We define it here default for human genome assembly(hg18,h19) to be 1e-6 as GATK did ", required = false)
    public double referenceGenomeErr = 1e-6;
	
	@Argument(fullName = "reference_genome_version", shortName = "rgv", doc = "Reference genome assembly version, the default value is hg18 ", required = false)
    public String referenceGenomeVer = "hg18";
	
	@Argument(fullName = "ti_vs_tv", shortName = "tvt", doc = "Transition rate vs. Transversion rate, in human genome, the default is 2", required = false)
    public int tiVsTv = 2;
	
	//@Argument(fullName = "tcga_format_vcf", shortName = "tcga", doc = "output TCGA specific VCF format or not, not used yet, in test", required = false)
    //public boolean tcga = false;
	
    @Argument(fullName = "output_modes", shortName = "out_modes", doc = "Output modes[EMIT_VARIANTS_ONLY,EMIT_ALL_CONFIDENT_SITES,EMIT_ALL_SITES,EMIT_ALL_CPG, EMIT_ALL_CYTOSINES,EMIT_HET_SNPS_ONLY, DEFAULT_FOR_TCGA]", required = false)
    public BisulfiteGenotyperEngine.OUTPUT_MODE OutputMode = BisulfiteGenotyperEngine.OUTPUT_MODE.DEFAULT_FOR_TCGA;
    
  //  @Argument(fullName = "vcf_file_name", shortName = "vfn", doc = "output Vcf file", required = true)
//	public String vfn = null;
    
    @Argument(fullName = "vcf_file_name_1", shortName = "vfn1", doc = "output Vcf file, when used for [DEFAULT_FOR_TCGA] output mode, it is used to store all CpG sites. While the original vcf file is to store all CpG sites", required = false)
	public String vfn1 = null;
    
    @Argument(fullName = "vcf_file_name_2", shortName = "vfn2", doc = "output Vcf file 2, only used for [DEFAULT_FOR_TCGA] output mode, it is used to store all SNP sites. While the original vcf file is to store all CpG sites", required = false)
	public String vfn2 = null;
	
    @Argument(fullName = "output_reads_after_downsampling", shortName = "orad", doc = "output Bam file that after downsapling, for performance test only", required = false)
    public boolean orad = false;
    
    @Argument(fullName = "file_name_output_reads_after_downsampling", shortName = "fnorad", doc = "output Bam file's name that after downsapling, for performance test only", required = false)
	public String fnorad = null;
    
    @Argument(fullName = "output_reads_coverage_after_downsampling", shortName = "orcad", doc = "define output Bam file's mean coverage that after downsapling, for performance test only", required = false)
	public int orcad = 1;
    
    @Argument(fullName = "file_name_output_cpg_reads_detail", shortName = "fnocrd", doc = "output CpG reads bed file that contain each CpG's position in reads information, for test only", required = false)
	public String fnocrd = null;
    
    @Argument(fullName = "file_name_output_verbose_detail", shortName = "fnovd", doc = "output file that contain verbose information, for test only", required = false)
	public String fnovd = null;
    
    @Argument(fullName = "output_verbose_detail", shortName = "ovd", doc = "output_verbose_detail, for performance test only", required = false)
    public boolean ovd = false;
    
    @Argument(fullName = "locus_not_continuous", shortName = "lnc", doc = "locu to look at is not continuous, if the distance is too large, it will make some trouble in multithread VCF writer, just enable this option in performance test only", required = false)
    public boolean lnc = false;
	
	public BisulfiteArgumentCollection clone() {
		BisulfiteArgumentCollection bac = new BisulfiteArgumentCollection();
		bac.GLmodel = GLmodel;
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
        bac.bisulfiteConversionModeOnestrand = bisulfiteConversionModeOnestrand;
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
        bac.minConv = minConv;
        bac.bsRate = bsRate;
        bac.overRate = overRate;
        bac.validateDbsnpHet = validateDbsnpHet;
        bac.novelDbsnpHet = novelDbsnpHet;
        bac.referenceGenomeErr = referenceGenomeErr;
        bac.referenceGenomeVer = referenceGenomeVer;
        bac.heterozygosity = heterozygosity;
        bac.tiVsTv = tiVsTv;
        
        bac.orad = orad;
        bac.fnorad = fnorad;
        bac.vfn2 = vfn2;
        bac.fnocrd = fnocrd;
        bac.vfn1 = vfn1;
        bac.fnovd = fnovd;
        bac.ovd = ovd;
        bac.lnc = lnc;
        
        return bac;
    }
	
	@Override
	public String toString() {
		String cls = new String();
		//int i = 0;
		for (Field f : this.getClass().getDeclaredFields()){
			cls = cls.concat("\n");
			cls = cls.concat(f.toString());
			//i++;
		}
		return cls;
		
		
	}

}
