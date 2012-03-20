package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;


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

public class BisulfiteVCFConstants{


	    // INFO/FORMAT field keys for Bis-SNP VCF 
		public static final String GENOTYPE_TYPE = "HOM_REF,HET,HOM_VAR";
	    public static final String NUMBER_OF_C_KEY = "Cm";
	    public static final String NUMBER_OF_T_KEY = "Cu";
	    public static final String C_STRAND_KEY = "Strand";
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
	    public static final String PROGRAM_ARGS = "BisSNP Program Args";
	    public static final String ID_KEY = "ID";

}

