package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import org.broad.tribble.vcf.VCFHeaderVersion;


public class BisulfiteVCFConstants{


	    // INFO/FORMAT field keys for Bis-SNP VCF 
		public static final String GENOTYPE_TYPE = "HOM_REF,HET,HOM_VAR";
	    public static final String NUMBER_OF_C_KEY = "NumC";
	    public static final String NUMBER_OF_T_KEY = "NumT";
	    public static final String C_IN_NEG_STRAND_KEY = "NegC";
	    public static final String CYTOSINE_TYPE = "CType";
	    public static final String CYTOSINE_METHY_VALUE = "methy";
	    public static final String VCF_HEADER_VERSION_FORMAT = "fileformat";
	    public static final String VCF_HEADER_VERSION_DATE = "fileDate";
	    public static final String VCF_HEADER_VERSION_TCGA_VERSION = "tcgaversion";
	    public static final String VCF_HEADER_VERSION_LOG = "vcfProcessLog";
	    public static final String VCF_HEADER_VERSION_REF = "reference";
	    public static final String VCF_HEADER_VERSION_ASSEMBLY = "assembly";
	    public static final String VCF_HEADER_VERSION_CENTER = "center";
	    public static final String VCF_HEADER_VERSION_PHASE = "phasing";
	    public static final String VCF_HEADER_VERSION_GAF = "geneAnno";
	    public static final String GENOTYPE_LIKELIHOODS_KEY = "GL";

}

