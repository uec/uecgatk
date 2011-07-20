package edu.usc.epigenome.uecgatk.bisulfitesnpmodel;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.ResourceBundle;

import net.sf.picard.filter.SamRecordFilter;

import org.broad.tribble.TribbleException;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.ArgumentTypeDescriptor;
import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.gatk.CommandLineExecutable;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.broadinstitute.sting.gatk.datasources.reads.SAMReaderID;
import org.broadinstitute.sting.gatk.io.stubs.OutputStreamArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileReaderArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet.RMDStorageType;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.gatk.walkers.Attribution;
import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.ApplicationDetails;
import org.broadinstitute.sting.utils.text.TextFormattingUtils;
import org.broadinstitute.sting.utils.text.XReadLines;

import org.broadinstitute.sting.gatk.uscec.bisulfitesnpmodel.BisulfiteGenotyper;

public class MethylSNP extends CommandLineExecutable {
	
	@Argument(fullName = "analysis_type", shortName = "T", doc = "Type of analysis to run")
    private String analysisName = null;
	
	@Argument(fullName = "auto_estimate_cytosine_mode", shortName = "aecm", doc = "the first run would be to run auto_estimate_cytosine methylation status")
    private static boolean autoEstimateC = false;

	//copy from GATK, since they are private class in GATK
	private final Collection<Object> bisulfiteArgumentSources = new ArrayList<Object>();
	
    // our argument collection, the collection of command line args we accept
    @ArgumentCollection
    private GATKArgumentCollection argCollection = new GATKArgumentCollection();
    
	//@ArgumentCollection
   // private BisulfiteArgumentCollection bisulfiteArgCollection = new BisulfiteArgumentCollection();
	
	private static boolean secondIteration = false;
	private static CytosineTypeStatus cts = null;
	
	public Walker<?,?> walker = null;
	
	protected VCFWriter writer = null;
	private VariantAnnotatorEngine annotationEngine = null;
	
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

   // protected BisulfiteArgumentCollection getBisulfiteArgumentCollection() {
  //      return bisulfiteArgCollection;
  // }

    

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			MethylSNP instance = new MethylSNP();
			
			
				start(instance, args);
				secondIteration = true;
				if(autoEstimateC & secondIteration){
					//instance.setupInfo();
					//start(instance, args);
					instance.execute();
					
					//start(instance, args);
				}

            System.exit(CommandLineProgram.result); // todo -- this is a painful hack

            
        } catch (UserException e) {
            exitSystemWithUserError(e);
        } catch (TribbleException e) {
            // We can generate Tribble Exceptions in weird places when e.g. VCF genotype fields are
            //   lazy loaded, so they aren't caught elsewhere and made into User Exceptions
            exitSystemWithUserError(e);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
        
	}
	
	
	public void setupInfo(){
		if(walker instanceof BisulfiteGenotyper){
			((BisulfiteGenotyper) walker).setWriter(writer);
			((BisulfiteGenotyper) walker).setAnnoEng(annotationEngine);
			//System.err.println("writer2: " + writer.toString());
			
		}
	}
	
	public static List<String> createApplicationHeader() {
        String version = "bs-0.15";
		List<String> header = new ArrayList<String>();
        header.add(String.format("The MethylSNP (BSSNP) v%s, Compiled %s",version, getBuildTime()));
        header.add(String.format("Based on The Genome Analysis Toolkit (GATK) v%s (in sorceforge tree, the version number is 5288)",getVersionNumber()));
        header.add("Copyright (c) 2011 USC Epigenome Center");
        header.add("Please view our documentation at http://wiki.epigenome.usc.edu/twiki/bin/view");
        header.add("For support, please send email to yapingli@usc.edu or benbfly@gmail.com");
        return header;
    }
	
	//need to figure out how to get version number
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
        	if(autoEstimateC & secondIteration){
        		System.out.println("2nd iteration!");
        		//engine = new GenomeAnalysisEngine();
        		
        		bisulfiteArgumentSources.clear();
        		//engine.setParser(parser);
        		bisulfiteArgumentSources.add(this);
        		
        		//engine.setArguments(getArgumentCollection());
        		
                // File lists can require a bit of additional expansion.  Set these explicitly by the engine. 
               // engine.setSAMFileIDs(unpackBAMFileList(getArgumentCollection()));
               // engine.setReferenceMetaDataFiles(unpackRODBindings(getArgumentCollection()));
                 walker = (BisulfiteGenotyper) engine.getWalkerByName(getAnalysisName());
        		((BisulfiteGenotyper) walker).setAutoParameters(autoEstimateC, secondIteration);
        		setupInfo();
        		engine.setWalker(walker);
                walker.setToolkit(engine);

              //  Collection<SamRecordFilter> filters = engine.createFilters();
               // engine.setFilters(filters);
        		((BisulfiteGenotyper) walker).setCytosineMethyStatus(cts);
        		//System.err.println("writer: " + writer.toString());
        		
                // load the arguments into the walker / filters.
                // TODO: The fact that this extra load call exists here when all the parsing happens at the engine
                // TODO: level indicates that we're doing something wrong.  Turn this around so that the GATK can drive
                // TODO: argument processing.
                loadArgumentsIntoObject(walker);
                bisulfiteArgumentSources.add(walker);
                Collection<SamRecordFilter> filters = engine.getFilters();
                for (SamRecordFilter filter: filters) {
                    loadArgumentsIntoObject(filter);
                    bisulfiteArgumentSources.add(filter);
                }
               // System.err.println(engine.getArguments().referenceFile.toString());
                
                engine.execute();
                
        	}
        	else{
        		//System.err.println("1st iteration!");
        		 engine.setParser(parser);
        	     bisulfiteArgumentSources.add(this);
        	        
        		walker = engine.getWalkerByName(getAnalysisName());
        		if(walker instanceof BisulfiteGenotyper){
        			
        			((BisulfiteGenotyper) walker).setAutoParameters(autoEstimateC, secondIteration);
        		}
        		
        		
        		engine.setArguments(getArgumentCollection());

                // File lists can require a bit of additional expansion.  Set these explicitly by the engine. 
                engine.setSAMFileIDs(unpackBAMFileList(getArgumentCollection()));
                engine.setReferenceMetaDataFiles(unpackRODBindings(getArgumentCollection()));

                engine.setWalker(walker);
                walker.setToolkit(engine);

                Collection<SamRecordFilter> filters = engine.createFilters();
                engine.setFilters(filters);
                

                // load the arguments into the walker / filters.
                // TODO: The fact that this extra load call exists here when all the parsing happens at the engine
                // TODO: level indicates that we're doing something wrong.  Turn this around so that the GATK can drive
                // TODO: argument processing.
                loadArgumentsIntoObject(walker);
                bisulfiteArgumentSources.add(walker);

                for (SamRecordFilter filter: filters) {
                    loadArgumentsIntoObject(filter);
                    bisulfiteArgumentSources.add(filter);
                }
                if(walker instanceof BisulfiteGenotyper){
        			writer = ((BisulfiteGenotyper) walker).getWriter();
        			//System.err.println("writer1: " + writer.toString());
        			this.annotationEngine = ((BisulfiteGenotyper) walker).getAnnoEng();
        		}
                
                engine.execute();
                if(walker instanceof BisulfiteGenotyper){
        			cts = ((BisulfiteGenotyper) walker).getCytosineMethyStatus();
        			for(String key : cts.cytosineListMap.keySet()){
        				Double[] values = cts.cytosineListMap.get(key);
        				for(Double value : values){
        					//System.err.println("cts.key: " + key + "\tcts.value: " + value);
        				}
        			}
        			
        			//System.err.println(cts.chhMethyLevel);
        			//System.err.println(cts.chgMethyLevel);
        			//System.err.println(cts.cpgMethyLevel);
        		}
                //System.err.println(result.toString());
                //generateGATKRunReport(walker);
                
        	}
 
        } catch ( Exception e ) {
           // generateGATKRunReport(walker, e);
            throw e;
        }

        // always return 0
        return 0;
    }

    
    //copy from GATK, since they are private class in GATK, or I could just limit the number of BAM used
    /**
     * Unpack the bam files to be processed, given a list of files.  That list of files can
     * itself contain entries which are lists of other files to be read (note: you cannot have lists of lists of lists)
     *
     * @param argCollection the command-line arguments from which to extract the BAM file list.
     * @return a flattened list of the bam files provided
     */
    private List<SAMReaderID> unpackBAMFileList(GATKArgumentCollection argCollection) {
        List<SAMReaderID> unpackedReads = new ArrayList<SAMReaderID>();
        for( String inputFileName: argCollection.samFiles ) {
            Tags inputFileNameTags = parser.getTags(inputFileName);
            inputFileName = expandFileName(inputFileName);
            if (inputFileName.toLowerCase().endsWith(".list") ) {
                try {
                    for(String fileName : new XReadLines(new File(inputFileName)))
                        unpackedReads.add(new SAMReaderID(fileName,parser.getTags(inputFileName)));
                }
                catch( FileNotFoundException ex ) {
                    throw new UserException.CouldNotReadInputFile(new File(inputFileName), "Unable to find file while unpacking reads", ex);
                }
            }
            else if(inputFileName.toLowerCase().endsWith(".bam")) {
                unpackedReads.add(new SAMReaderID(inputFileName,inputFileNameTags));
            }
            else if(inputFileName.endsWith("stdin")) {
                unpackedReads.add(new SAMReaderID(inputFileName,inputFileNameTags));
            }
            else {
                throw new UserException.CommandLineException(String.format("The GATK reads argument (-I) supports only BAM files with the .bam extension and lists of BAM files " +
                        "with the .list extension, but the file %s has neither extension.  Please ensure that your BAM file or list " +
                        "of BAM files is in the correct format, update the extension, and try again.",inputFileName));
            }
        }
        return unpackedReads;
    }
    
  //copy from GATK, since they are private class in GATK
    /**
     * Convert command-line argument representation of ROD bindings to something more easily understandable by the engine.
     * @param argCollection input arguments to the GATK.
     * @return a list of expanded, bound RODs.
     */
    private Collection<RMDTriplet> unpackRODBindings(GATKArgumentCollection argCollection) {
        Collection<RMDTriplet> rodBindings = new ArrayList<RMDTriplet>();

        for (String fileName: argCollection.RODBindings) {
            Tags tags = parser.getTags(fileName);
            fileName = expandFileName(fileName);

            List<String> positionalTags = tags.getPositionalTags();
            if(positionalTags.size() != 2)
                throw new UserException("Invalid syntax for -B (reference-ordered data) input flag.  " +
                        "Please use the following syntax when providing reference-ordered " +
                        "data: -B:<name>,<type> <filename>.");
            // Assume that if tags are present, those tags are name and type.
            // Name is always first, followed by type.
            String name = positionalTags.get(0);
            String type = positionalTags.get(1);

            RMDStorageType storageType = null;
            if(tags.getValue("storage") != null)
                storageType = Enum.valueOf(RMDStorageType.class,tags.getValue("storage"));
            else if(fileName.toLowerCase().endsWith("stdin"))
                storageType = RMDStorageType.STREAM;
            else
                storageType = RMDStorageType.FILE;

            rodBindings.add(new RMDTriplet(name,type,fileName,storageType));
        }

        if (argCollection.DBSNPFile != null) {
            if(argCollection.DBSNPFile.toLowerCase().contains("vcf"))
                throw new UserException("--DBSNP (-D) argument currently does not support VCF.  To use dbSNP in VCF format, please use -B:dbsnp,vcf <filename>.");

            String fileName = expandFileName(argCollection.DBSNPFile);
            RMDStorageType storageType = fileName.toLowerCase().endsWith("stdin") ? RMDStorageType.STREAM : RMDStorageType.FILE;

            rodBindings.add(new RMDTriplet(DbSNPHelper.STANDARD_DBSNP_TRACK_NAME,"dbsnp",fileName,storageType));
        }

        return rodBindings;
    }
    
    private String expandFileName(String argument) {
        if(argument.trim().equals("-"))
            return "/dev/stdin";
        return argument;
    }
    
  //copy from GATK, since they are private class in GATK
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
    
    
    
}
