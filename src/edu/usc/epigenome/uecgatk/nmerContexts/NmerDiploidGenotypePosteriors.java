package edu.usc.epigenome.uecgatk.nmerContexts;

/**
 * @author benb
 * 
 * This is output by BisSNP or other genotype callers to give
 * posterior probabilities for each position in the context.
 * It may be a good idea to make this class contain multiple
 * read groups, because then we could put functions in here
 * to determine whether there is a single context across read
 * groups or whether there are several.  Those functions could
 * also go in the class managing these objects (for instance
 * BisulfitePileup)
 * 
 *
 */
public class NmerDiploidGenotypePosteriors {

}
