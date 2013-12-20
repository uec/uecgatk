package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Map;

import org.broadinstitute.sting.commandline.Argument;

import edu.usc.epigenome.uecgatk.benWalkers.CpgWalker;
import edu.usc.epigenome.uecgatk.benWalkers.CpgWalkerWindowWigWriter;
import edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.BisulfiteSeqToCytosineWindowWalker.ContextConditions;

public class BisulfiteSeqToCytosineVariableWigWalker extends
		BisulfiteSeqToCytosineWindowWalker {

    @Argument(fullName = "csvMode", shortName = "csv", doc = "Output csv rather than bedgraph", required = false)
    public boolean csvMode = false;

	/**
	 * @param cpgWalkerType
	 */
	public BisulfiteSeqToCytosineVariableWigWalker() {
		super("edu.usc.epigenome.genomeLibs.MethylDb.CpgWalker.CpgWalkerWindowWigWriter");
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.BisulfiteSeqToCytosineWindowWalker#initialize()
	 */
	@Override
	public void initialize() {
		// TODO Auto-generated method stub
		super.initialize();
		
		for (ContextConditions cond : this.validContexts())
		{
			CpgWalkerWindowWigWriter walker = (CpgWalkerWindowWigWriter)this.walkerByCondition.get(cond);
			
			String maxsec = (this.useFixedWind) ? String.format("fixedWind%d", this.fixedWindMinLength) : String.format("varMaxStretch%d", this.maxWindStretch);
			
			String name = String.format("%s-%s-minc%d-minct%d-%s", this.outPrefix, cond.name(), this.minCpgs, this.minCT, maxsec);
			String fn = String.format("%s.wig", name);
			String fnBedgraph = String.format("%s.%s", name, (csvMode)?"csv":"bedGraph");
			System.err.printf("Initializing file: %s\n",fn);
			PrintWriter pw = null;
			PrintWriter pwBedgraph = null;
			try
			{
				pw = new PrintWriter(new FileOutputStream(fn));
				pwBedgraph = new PrintWriter(new FileOutputStream(fnBedgraph));
			}
			catch (Exception e)
			{
				System.err.printf("Fatal error, CpgWalkerWindowWigWriter couldn't write to file \"%s\"\n%s",fn,e.toString());
				e.printStackTrace();
				System.exit(1);
			}
			String desc = String.format("%s-wig", name);
			pw.printf("track type=wiggle_0 name=%s description=%s color=204,102,0 visibility=full " +
					" graphType=points autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=0:100\n", desc, desc);	
			walker.setWigWriter(pw);

			if (!csvMode)
			{
				desc = String.format("%s-bg", name);
				pwBedgraph.printf("track type=bedGraph name=%s description=%s color=204,102,0 visibility=full " +
						" graphType=points autoScale=off alwaysZero=off maxHeightPixels=64:32:10 viewLimits=-1:100\n", desc, desc);
			}
			walker.setBedgraphWriter(pwBedgraph);
			if (csvMode) walker.setCsvMode(true);
		
		}
	}

	/* (non-Javadoc)
	 * @see edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers.BisulfiteSeqToCytosineWindowWalker#onTraversalDone(java.util.Map)
	 */
	@Override
	public void onTraversalDone(Map<ContextConditions, CpgWalker> result) {
		super.onTraversalDone(result);

		for (ContextConditions cond : this.validContexts())
		{
			CpgWalkerWindowWigWriter walker = (CpgWalkerWindowWigWriter)this.walkerByCondition.get(cond);
			walker.getWigWriter().close();
			walker.getBedgraphWriter().close();
		}
	
	}
	
	
	

}
