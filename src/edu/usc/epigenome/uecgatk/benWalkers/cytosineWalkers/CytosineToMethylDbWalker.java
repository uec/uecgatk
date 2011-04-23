package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.File;
import java.io.FileOutputStream;
import java.nio.CharBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.nio.charset.CharsetEncoder;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.locks.ReentrantLock;


import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;



public class CytosineToMethylDbWalker extends LocusWalkerToBisulfiteCytosineWalker<Integer, Long> {


    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////

    @Argument(fullName = "outPrefix", shortName = "pre", doc = "Output prefix for all output files", required = true)
    public String outPrefix = null;
	
    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
	protected static Map<String,FileChannel> outfilesByChr = new HashMap<String,FileChannel>();
	private static final ReentrantLock outfilesMapLock = new ReentrantLock();
	protected Charset charset = Charset.forName("UTF-8");
	//protected CharsetEncoder encoder = charset.newEncoder();  // If the encoder is static we have concurrency issues.
	
    /////////////////////////////
    // Static Variables
    /////////////////////////////


	//---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

	
	//---------------------------------------------------------------------------------------------------------------
    //
    // LocusWalker overrides
    //
    //---------------------------------------------------------------------------------------------------------------

	/**
	 * Provides an initial value for the reduce function.  Hello walker counts loci,
	 * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
	 * @return 0.
	 */
	@Override
	public Long reduceInit()
	{ 
		Long out = 0L;
		return out;
	}



	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 * 
	 * Must override this to get multithreaded capability
	 */
	@Override
	public Long treeReduce(Long lhs, Long rhs) 
	{
		return lhs + rhs;
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Long result) 
	{
		out.printf("Saw %d cytosines\n", result);
		for (String chr : outfilesByChr.keySet())
		{
			FileChannel fc = outfilesByChr.get(chr);
			try
			{
				fc.close();
			}
			catch (Exception e)
			{
				System.err.printf("Fatal error , could close file for %s\n%s\n", chr, e.toString());
				e.printStackTrace();
				System.exit(1);
			}
		}
	}

	
	/***************************************************
	 * cytosine walker overrides
	 ***************************************************/
	
	
	@Override
	protected void alertNewContig(String newContig) 
	{
	}

	
	@Override
	protected Integer processCytosine(CpgBackedByGatk thisC)
	{
		// Get file
		
		FileChannel fc = null;
		String chr = thisC.getChrom();
		
		// Threads can race and create multiple files.
		outfilesMapLock.lock();
		try{
			if (outfilesByChr.containsKey(chr))
			{
				fc = outfilesByChr.get(chr);
				//logger.info("Found wig " + wig.toString());
			}
			else
			{
				String name = String.format("%s-%s", this.outPrefix, chr);
				String outfn = String.format("%s.txt",name);

				logger.info("NEW file " + outfn);
					FileOutputStream fout=null;

					fout = new FileOutputStream(outfn);
					fc = fout.getChannel();
					outfilesByChr.put(chr, fc);
			}
		}
		catch (Exception e)
		{
			System.err.printf("Fatal error , could not write to file for  %s\n%s\n", chr, e.toString());
			e.printStackTrace();
			System.exit(1);
		}
		finally
		{
			outfilesMapLock.unlock();
		}
		
		


		/**
		 * @author benb
		 * 
		 * Generates MethylDb table:
		 * 
		   CREATE TABLE `TEMPLATE_CHR` (
		  `chromPos` INT UNSIGNED NOT NULL,
		  `strand` enum('+','-') NOT NULL,
		  `totalReads` SMALLINT UNSIGNED NOT NULL,
		  `cReads` SMALLINT UNSIGNED NOT NULL,
		  `cReadsNonconversionFilt` SMALLINT UNSIGNED NOT NULL,
		  `tReads` SMALLINT UNSIGNED NOT NULL,
		  `agReads` SMALLINT UNSIGNED NOT NULL,

		  `totalReadsOpposite` SMALLINT UNSIGNED NOT NULL,
		  `aReadsOpposite` SMALLINT UNSIGNED NOT NULL,

		  `nextBaseGreads` SMALLINT UNSIGNED NULL default '0',
		  `nextBaseTotalReads` SMALLINT UNSIGNED NULL default '0',

		  `cpgWeight` SMALLINT UNSIGNED NOT NULL,

		  -- Must have a primary key to use updatable rows
		  PRIMARY KEY chromPos(chromPos)
		 )
			
		 * 
		 */

		// Get cytosine info
//		String context = thisC.context();
//		AlignmentContext ac = thisC.getAlignmentContext();
//		double meth = thisC.fracMeth(false);
//		int pos = thisC.chromPos;
		
		try
		{
			
			
			fc.write(charset.newEncoder().encode(CharBuffer.wrap( thisC.toString() + "\n" )));
		}
		catch (Exception e)
		{
			System.err.printf("Fatal error , could not write to file for %s\n%s\n", chr, e.toString());
			e.printStackTrace();
			System.exit(1);
		}
		

		return 1;
	}





	@Override
	protected Long reduceCytosines(Integer value, Long sum) 
	{
		return value.longValue() + sum;
	}
	
}
