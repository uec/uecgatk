package edu.usc.epigenome.uecgatk.BisSNP;

import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.broadinstitute.sting.commandline.Advanced;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.commandline.RodBinding;

import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyperEngine;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteEnums.cytosineContextParameters;

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
	
	 @Input(fullName="dbsnp", shortName = "D", doc="dbSNP file", required=false)
	 public RodBinding<VariantContext> dbsnp;
	
	@Argument(fullName = "sequencing_mode", shortName = "sm", doc = "Bisulfite mode: BM, GNOMe-seq mode: GM, Normal sequencing mode: NM", required = false)
    public BisulfiteEnums.MethylSNPModel sequencingMode = BisulfiteEnums.MethylSNPModel.BM;
	
//	@Argument(fullName = "paired_end_mode", shortName = "pem", doc = "work in paired end mode", required = false)
  //  public boolean pairedEndMode = false;
	
	@Argument(fullName = "bisulfite_conversion_only_on_one_strand", shortName = "bcm", doc = "true: Illumina protocol which is often used, only bisulfite conversion strand is kept (Lister protocol, sequence 2 forward strands only); false: Cokus protocol, sequence all 4 bisulfite converted strands", required = false)
    public boolean bisulfiteConversionModeOnestrand = true;
	
	@Advanced
	@Argument(fullName = "cytosine_contexts_acquired", shortName = "C", doc = "Specify the cytosine contexts to check (e.g. -C CG,1,0.7  CG is the methylation pattern to check, 1 is the C's position in CG pattern, 0.7 is the forced cytosine pattern methylation pattern. You could specify '-C' multiple times for different cytosine pattern)", required = false)
    public List<String> cytosineContextsAcquired = new ArrayList<String>();

	@Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum mapping quality required to consider a base for calling", required = false)
    public int MIN_MAPPING_QUALTY_SCORE = 30;
	
	@Argument(fullName = "use_badly_mated_reads", shortName = "badMate", doc = "use badly mated reads for calling", required = false)
    public boolean USE_BADLY_MATED_READS = false;
	
	@Argument(fullName = "assume_single_sample", shortName = "single_sample", doc = "only single sample", required = false)
    public String ASSUME_SINGLE_SAMPLE = null;
	
	 @Argument(fullName = "max_mismatches", shortName = "mm40", doc = "Maximum number of mismatches within a 40 bp window (20bp on either side) around the target position for a read to be used for calling", required = false)
	 public int MAX_MISMATCHES = 3;

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
    public BisulfiteEnums.OUTPUT_MODE OutputMode = BisulfiteEnums.OUTPUT_MODE.DEFAULT_FOR_TCGA;
    
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
    
    @Argument(fullName = "bissnp_methy_summary_file", shortName = "bmsf", doc = "input the methylation summary estimate from BisSNP, for BisSNPUtils right now only", required = false)
	public String bmsf = null;
    
    public CytosinePatternsUserDefined cytosineDefined = null;
    
   // public BisulfiteArgumentCollection(){
    	//System.err.println(cytosineContextsAcquired.toString() + "\t" + cytosineContextsAcquired.isEmpty());
    	//System.err.println(coverageThresholds[0] + coverageThresholds[1] );
    	//cytosineDefined = new CytosinePatternsUserDefined(cytosineContextsAcquired);
    //	makeCytosine();
    	
   // }
    
    public void makeCytosine(){
    	cytosineDefined = new CytosinePatternsUserDefined(cytosineContextsAcquired);
  //  	System.err.println(cytosineContextsAcquired.toString() + "\t" + cytosineDefined.getContextDefined().keySet().toString());
    }
	
	public BisulfiteArgumentCollection clone() {
		BisulfiteArgumentCollection bac = new BisulfiteArgumentCollection();
		
		bac.GLmodel = GLmodel;
        bac.PCR_error = PCR_error;
        bac.GenotypingMode = GenotypingMode;
        bac.OutputMode = OutputMode;
        bac.NO_SLOD = NO_SLOD;
        //bac.ASSUME_SINGLE_SAMPLE = ASSUME_SINGLE_SAMPLE;
        bac.STANDARD_CONFIDENCE_FOR_CALLING = STANDARD_CONFIDENCE_FOR_CALLING;
        bac.STANDARD_CONFIDENCE_FOR_EMITTING = STANDARD_CONFIDENCE_FOR_EMITTING;
        bac.MIN_BASE_QUALTY_SCORE = MIN_BASE_QUALTY_SCORE;
        bac.MIN_MAPPING_QUALTY_SCORE = MIN_MAPPING_QUALTY_SCORE;
        bac.MAX_MISMATCHES = MAX_MISMATCHES;
        bac.USE_BADLY_MATED_READS = USE_BADLY_MATED_READS;
        bac.MAX_DELETION_FRACTION = MAX_DELETION_FRACTION;
        bac.MIN_INDEL_COUNT_FOR_GENOTYPING = MIN_INDEL_COUNT_FOR_GENOTYPING;
        bac.INDEL_HETEROZYGOSITY = INDEL_HETEROZYGOSITY;
        
        bac.dbsnp = dbsnp;
        bac.sequencingMode = sequencingMode;
 //       bac.pairedEndMode = pairedEndMode;
        bac.bisulfiteConversionModeOnestrand = bisulfiteConversionModeOnestrand;
        bac.cytosineContextsAcquired = cytosineContextsAcquired;
        //bac.cytosineDefined = cytosineDefined;
        bac.makeCytosine();
        //bac.cytosineDefined = cytosineDefined;
        
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
        bac.bmsf = bmsf;
        
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
