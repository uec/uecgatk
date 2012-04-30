package edu.usc.epigenome.uecgatk.BisSNP;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.ResourceBundle;

import net.sf.picard.filter.SamRecordFilter;


import org.broadinstitute.sting.utils.codecs.vcf.SortingVCFWriter;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.arguments.DbsnpArgumentCollection;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.filters.ReadFilter;
import org.broadinstitute.sting.gatk.io.stubs.OutputStreamArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileReaderArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet.RMDStorageType;
import org.broadinstitute.sting.gatk.walkers.Attribution;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.ApplicationDetails;
import org.broadinstitute.sting.utils.text.ListFileUtils;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.RodBinding;

import edu.usc.epigenome.uecgatk.BisSNP.BisulfiteGenotyper;

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

public class BisSNP extends CommandLineExecutable {
	
	@Argument(fullName = "analysis_type", shortName = "T", doc = "Type of analysis to run")
    private String analysisName = null;
	
	@Argument(fullName = "auto_estimate_cytosine_mode", shortName = "aecm", doc = "automately estimate cytosine pattern methylation status in the first iteration")
    private static boolean autoEstimateC = false;
		
	 //control the output, output to TCGA VCF 
    //@Output(doc="File to which variants should be written",required=false)
   // protected SortingVCFWriter writer = null;
 
	private static String BisVersion = "BisSNP-0.63";
	
	private final Collection<Object> bisulfiteArgumentSources = new ArrayList<Object>();
	
    // argument collection, the collection of command line args we accept. copy from GATK, since they are private class in GATK
    @ArgumentCollection
    private GATKArgumentCollection argCollection = new GATKArgumentCollection();
    
    
    private static String argCommandline = "";
    
	//to record it is in second iteration or not
	private static boolean secondIteration = false;
	
	//to record cytosine pattern methylation status estimated in the first iteration
	private static CytosinePatternsUserDefined cytosineDefinedMemorizedForSecondRun = null;
	
	//private RodBinding<VariantContext> dbsnp = null;
	
	public Walker<?,?> walker = null;
		
	
	@Override
    protected ApplicationDetails getApplicationDetails() {
        return new ApplicationDetails(createApplicationHeader(),
        		Collections.<String>emptyList(),
                ApplicationDetails.createDefaultRunningInstructions(getClass()),
                null);
    }

	@Override
	public String getAnalysisName() {
		// TODO Auto-generated method stub
		return analysisName;
	}

	

	@Override
	protected GATKArgumentCollection getArgumentCollection() {
		// TODO Auto-generated method stub
		return argCollection;
	}

     	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			BisSNP instance = new BisSNP();
	        for (String str : args) {
	        	argCommandline = argCommandline + str + " ";
	        }
			//System.err.println(output);
			start(instance, args);
			secondIteration = true;
			if(autoEstimateC & secondIteration){ 
				instance.execute(); // do the second iteration if it is in two-iteration mode
			}
            System.exit(CommandLineProgram.result);      
        } catch (UserException e) {
            exitSystemWithUserError(e);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
        
	}
	
	//set up Writer information. if writer is not initiat here, then there will be some wired close stream problem.
	public void setupInfo(){
		if(walker instanceof BisulfiteGenotyper){
			
			if(argCollection.numberOfThreads == 1){
			//	((BisulfiteGenotyper) walker).setWriter(writer);

			}
			else{
			//	((BisulfiteGenotyper) walker).setWriter(writer);
			}
		}
	}
	
	public static List<String> createApplicationHeader() {

		List<String> header = new ArrayList<String>();
        header.add(String.format("The %s, Compiled %s", getBisSNPVersionNumber(), getBuildTime()));
        header.add(String.format("Based on The Genome Analysis Toolkit (GATK) v%s (prebuild GATK package could be download here: ftp://ftp.broadinstitute.org/pub/gsa/GenomeAnalysisTK/GenomeAnalysisTK-1.5-3-gbb2c10b.tar.bz2)",getVersionNumber()));
        header.add("Copyright (c) 2011 USC Epigenome Center");
        header.add("Please view our documentation at http://epigenome.usc.edu/publicationdata/bissnp2011/");
        header.add("For support, please send email to lyping1986@gmail.com or benbfly@gmail.com");
        return header;
    }

	public static String getVersionNumber() {
        ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
        
        return headerInfo.containsKey("org.broadinstitute.sting.gatk.version") ? headerInfo.getString("org.broadinstitute.sting.gatk.version") : "<unknown>";
    }
	
 
    public static String getBuildTime() {
        ResourceBundle headerInfo = TextFormattingUtils.loadResourceBundle("StingText");
        return headerInfo.containsKey("build.timestamp") ? headerInfo.getString("build.timestamp") : "<unknown>";
    }

    @Override
    protected int execute() throws Exception {
 
        try {
        	if(autoEstimateC & secondIteration){ // if auto-estimate cytosine model, and second iteration
        	//	System.out.println("2nd iteration!");

        		bisulfiteArgumentSources.clear();
        		
        		bisulfiteArgumentSources.add(this);
        		
                 walker = (BisulfiteGenotyper) engine.getWalkerByName(getAnalysisName());
        		((BisulfiteGenotyper) walker).setAutoParameters(autoEstimateC, secondIteration, argCommandline);
        	//	((BisulfiteGenotyper) walker).setBoundDbSNP(dbsnp);
        		setupInfo();
        		engine.setWalker(walker);
                walker.setToolkit(engine);

        		((BisulfiteGenotyper) walker).setCytosineMethyStatus(cytosineDefinedMemorizedForSecondRun); //transfer cytosine pattern methylation status estimated in the first iteration to 2nd iteration
        		
                loadArgumentsIntoObject(walker);
                bisulfiteArgumentSources.add(walker);
                
                
                
                Collection<ReadFilter> filters = engine.getFilters();
                for (SamRecordFilter filter: filters) {
                    loadArgumentsIntoObject(filter);
                    bisulfiteArgumentSources.add(filter);
                }
                
                
                engine.execute();
                
        	}
        	else{ //1st iteration
        		
        		 engine.setParser(parser);
        	     bisulfiteArgumentSources.add(this);
        	        
        		walker = engine.getWalkerByName(getAnalysisName());
        		if(walker instanceof BisulfiteGenotyper){
        			
        			((BisulfiteGenotyper) walker).setAutoParameters(autoEstimateC, secondIteration, argCommandline);
        		}
	
        		engine.setArguments(getArgumentCollection());
 
                engine.setSAMFileIDs(ListFileUtils.unpackBAMFileList(getArgumentCollection().samFiles,parser));
               // engine.setReferenceMetaDataFiles(unpackRODBindings(getArgumentCollection()));

                engine.setWalker(walker);
                walker.setToolkit(engine);

                Collection<ReadFilter> filters = engine.createFilters();
                engine.setFilters(filters);
                
                loadArgumentsIntoObject(walker);
                bisulfiteArgumentSources.add(walker);
                
                

                Collection<RMDTriplet> rodBindings = ListFileUtils.unpackRODBindings(parser.getRodBindings(), parser);

                // todo: remove me when the old style system is removed
                if ( getArgumentCollection().RODBindings.size() > 0 ) {
                    logger.warn("################################################################################");
                    logger.warn("################################################################################");
                    logger.warn("Deprecated -B rod binding syntax detected.  This syntax has been eliminated in GATK 1.2.");
                    logger.warn("Please use arguments defined by each specific walker instead.");
                    for ( String oldStyleRodBinding : getArgumentCollection().RODBindings ) {
                        logger.warn("  -B rod binding with value " + oldStyleRodBinding + " tags: " + parser.getTags(oldStyleRodBinding).getPositionalTags());
                    }
                    logger.warn("################################################################################");
                    logger.warn("################################################################################");
                    System.exit(1);
                }

                engine.setReferenceMetaDataFiles(rodBindings);

                for (SamRecordFilter filter: filters) {
                    loadArgumentsIntoObject(filter);
                    bisulfiteArgumentSources.add(filter);
                }
                setupInfo();
                //if(walker instanceof BisulfiteGenotyper){      			
        				//((BisulfiteGenotyper) walker).setWriter(writer);

        		//}
                
                engine.execute();
                if(walker instanceof BisulfiteGenotyper){
                	cytosineDefinedMemorizedForSecondRun = ((BisulfiteGenotyper) walker).getCytosineMethyStatus(); //receive cytosine pattern methylation status estimated in the first iteration
             //   	dbsnp = ((BisulfiteGenotyper) walker).getBoundDbSNP();
                	
        		}
                
                
        	}
 
        } catch ( Exception e ) {
            throw e;
        }
        return 0;
    }

    
    
    //copy from GATK, since they are private class in GATK 
    private String expandFileName(String argument) {
        if(argument.trim().equals("-"))
            return "/dev/stdin";
        return argument;
    }
    

    /**
     * Subclasses of CommandLinePrograms can provide their own types of command-line arguments.
     * @return A collection of type descriptors generating implementation-dependent placeholders.
     */
    @Override
    protected Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
        return Arrays.asList( new VCFWriterArgumentTypeDescriptor(engine,System.out,bisulfiteArgumentSources),
                              new SAMFileReaderArgumentTypeDescriptor(engine),
                              new SAMFileWriterArgumentTypeDescriptor(engine,System.out),
                              new OutputStreamArgumentTypeDescriptor(engine,System.out) );
    }
    
    public static String getBisSNPVersionNumber(){
    	return BisVersion;
    }
    
   
}
