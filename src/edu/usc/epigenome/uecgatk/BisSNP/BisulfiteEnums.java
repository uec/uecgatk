/**
 * 
 */
package edu.usc.epigenome.uecgatk.BisSNP;

/**
 * @author yaping
 * @contact lyping1986@gmail.com
 * @time 2012 Mar 19, 2012 4:09:03 PM
 * 
 */
public class BisulfiteEnums {
	
	
	public enum cytosineContextParameters{
		cytosinePosition,
		cytosineMethylation,
		cytosineStrand,
		methylationAutoestimated
	}
	
	public enum MethylSNPModel {
        BM,
        GM,
        NM
    }
	
	public enum OUTPUT_MODE {
        EMIT_VARIANTS_ONLY, //only confident variants
        EMIT_ALL_CONFIDENT_SITES,
        EMIT_ALL_SITES,
        EMIT_ALL_CPG, //only confident cpgs
        EMIT_ALL_CYTOSINES, //only confident cytosines
        EMIT_HET_SNPS_ONLY, //only confident heterozygous snps
        DEFAULT_FOR_TCGA, //output two vcf files: 1. confident variants vcf file (not contain cytosine info), 2. confident cpgs vcf file (only contain homozygous cpg info), 3. in future, output cpg_read files
        EMIT_VARIANT_AND_CYTOSINES //output all cytosines, and all variants, for performance test usage
    }
	
	
}
