package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.FileOutputStream;
import java.io.PrintWriter;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.BaseUtils;
import org.usckeck.genome.ChromFeatures;

import edu.usc.epigenome.genomeLibs.MatUtils;
import edu.usc.epigenome.uecgatk.BaseUtilsMore;
import edu.usc.epigenome.uecgatk.benWalkers.CpgBackedByGatk;
import edu.usc.epigenome.uecgatk.benWalkers.LocusWalkerToBisulfiteCytosineWalker;


/**
 * @author benb
 * 
 * This is pretty efficient in multi-threaded mode.  A test run took 2:00 minutes in -nt 1 mode and 0:30 in -nt 10 mode.
 * I worry about memory and Garbage Collection though, I haven't tested it thoroughly.
 * 
 */
public class MethylatedSequenceMatrixCytosineWalker extends LocusWalkerToBisulfiteCytosineWalker<double[][], double[][]> {

    @Argument(fullName = "matOutfilePrefix", shortName = "outpre", doc = "Prefix for the enoLOGOs matrix file", required = true)
    public String matOutfilePrefix = null;

    @Argument(fullName = "cContext", shortName = "context", doc = "context to pull (default HCH)", required = true)
    public String cContext = "HCH";

    @Argument(fullName = "outputFasta", shortName = "fasta", doc = "Fasta output to stdout", required = false)
    public boolean outputFasta = false;

    @Argument(fullName = "flankBasesUp", shortName = "up", doc = "Number of flanking bases upstream (default 5)", required = false)
    public Integer flankBasesUp = 5;

    @Argument(fullName = "flankBasesDown", shortName = "down", doc = "Number of flanking bases downstream (default 5)", required = false)
    public Integer flankBasesDown = 5;
  
    @Argument(fullName = "minMethReads", shortName = "minmeth", doc = "Number of reads with a C (default 2)", required = false)
    public Integer minMethReads = 2;
    
    /**
	 * locus walker overrides
	 */

    @Override
	public void initialize() {
		// TODO Auto-generated method stub
		super.initialize();
			
		//logger.info(String.format("Initializing map (%d)",walkerByCondition.hashCode()));
		this.outputCph = true; // Because GNOME-seq used GCH
	}       

	/**
	 * Provides an initial value for the reduce function.  Hello walker counts loci,
	 * so the base case for the inductive step is 0, indicating that the walker has seen 0 loci.
	 * @return 0.
	 */
	@Override
	public double[][] reduceInit()
	{ 
		double[][] out = new double[BaseUtils.BASES.length][this.flankBasesDown+this.flankBasesUp+1];
		return out;
	}



	/* (non-Javadoc)
	 * @see org.broadinstitute.sting.gatk.walkers.TreeReducible#treeReduce(java.lang.Object, java.lang.Object)
	 * 
	 * Must override this to get multithreaded capability
	 */
	@Override
	public double[][] treeReduce(double[][] a, double[][] b) 
	{
		if (this.outputFasta)
		{
			System.err.println("Use -nt 1 with --outputFasta");
			System.exit(1);
		}
		
		double[][] out = MatUtils.sumMats(a, b);
		return out;
	}

	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(double[][] result) 
	{

		if (this.matOutfilePrefix != null)
		{
			String outfn = String.format("ENOLOGO-%s-%s-%dminmeth.txt", this.matOutfilePrefix, this.cContext, this.minMethReads);
			try
			{
				PrintWriter pw = new PrintWriter(new FileOutputStream(outfn));
				for (int i = -1; i < result.length; i++)
				{
					for (int j = -1; j < result[0].length; j++)
					{
						if ((i==-1) && (j==-1))
						{
							pw.append("PO");
						}
						else if (i==-1)
						{
							pw.append(String.format("%d", j+1));
						}
						else if (j==-1)
						{
							pw.append((char)BaseUtils.baseIndexToSimpleBase(i));
						}
						else
						{
							pw.append(String.format("%d",(int)result[i][j]));
						}
						
						if (j<(result[0].length)) pw.append("\t");
					}
					pw.append("\n");
				}
				pw.close();
			}
      		catch (Exception e)
       		{
				System.err.printf("Could not write to LOGO file %s:\n%s\n",
						outfn,e.toString());
				e.printStackTrace();
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
	protected double[][] processCytosine(CpgBackedByGatk thisC)
	{
		String context = thisC.context();
		double[][] outMat = this.reduceInit();
		if (context.equalsIgnoreCase(this.cContext))
		{
			int nmeth = thisC.cReads;
			//double meth = thisC.fracMeth(false);
			if (nmeth >= this.minMethReads)
			{
				boolean revStrand = thisC.negStrand;
				int cpos = thisC.chromPos;
				int start = cpos - (revStrand ? this.flankBasesDown : this.flankBasesUp);
				int end = cpos + (revStrand ? this.flankBasesUp : this.flankBasesDown);
				String refBasesStr = null;
	    		try
	    		{
	    			byte[] refBases = BaseUtilsMore.toUpperCase(this.getToolkit().getReferenceDataSource().getReference().getSubsequenceAt(this.prevContig, start, end).getBases());
	    			if (revStrand) refBases = BaseUtils.simpleReverseComplement(refBases);
	    			StringBuffer buff = new StringBuffer(refBases.length);
	    			for (int j = 0; j < refBases.length; j++)
	    			{
	    				buff.append((char)refBases[j]);
	    				
	    				// And increment mat
	    				outMat[BaseUtils.simpleBaseToBaseIndex(refBases[j])][j]++;
	    			}
	    			refBasesStr = buff.toString();
	    		}
	      		catch (Exception e)
	       		{
					System.err.printf("Non-fatal error, could not getToolkit()getReferenceDataSource().getReference().getSubsequenceAt(%s,%d,%d)\n%s\n",
							this.prevContig, start, end,e.toString());
					e.printStackTrace();
	       		}
	      		
	      		if (this.outputFasta && (refBasesStr != null))
	      		{
	      			String name = String.format("> %s-%d%s:%dm", this.prevContig, cpos, (revStrand)?"-":"+",nmeth);
	      			out.append(name);
	      			out.append("\n");
	      			out.append(refBasesStr);
	      			out.append("\n");
	      		}
			}
			
		}
		
		return outMat;
	}





	@Override
	protected double[][] reduceCytosines(double[][] value, double[][] sum) 
	{
		double[][] out = MatUtils.sumMats(value, sum);
		return out;
	}
	
}
