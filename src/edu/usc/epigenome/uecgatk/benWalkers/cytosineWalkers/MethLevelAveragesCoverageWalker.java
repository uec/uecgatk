package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.ReadMetrics;

import edu.usc.epigenome.genomeLibs.ListUtils;
import edu.usc.epigenome.genomeLibs.MethylDb.CpgSummarizers.CpgMethLevelSummarizer;

public class MethLevelAveragesCoverageWalker extends MethLevelAveragesWalker {

    @Argument(fullName = "header", shortName = "header", doc = "Include header", required = false)
    public boolean includeHeader = false;

    public MethLevelAveragesCoverageWalker() {
		// TODO Auto-generated constructor stub
	}

	
	/**
	 * Retrieves the final result of the traversal.
	 * @param result The ultimate value of the traversal, produced when map[n] is combined with reduce[n-1]
	 *               by the reduce function. 
	 */
	@Override
	public void onTraversalDone(Map<String, CpgMethLevelSummarizer> result) 
	{
		List<String> header = new ArrayList();				
		List<String> list = new ArrayList();				
		
		// The built in functions for read filters didn't work so we have to get them this way.
		ReadMetrics rm = this.getToolkit().getCumulativeMetrics();
		long totalReads =  rm.getNumReadsSeen();
		Map<Class,Long> counts = rm.getCountsByFilter();
		for (Class c : counts.keySet())
		{
			long numbad = counts.get(c);
			totalReads -= numbad;
			//out.printf("\t%s: %d-%d=%d\n", c, totalReads+numbad, numbad,totalReads);
		}
		
		header.add("NumReads");
		list.add(String.format("%d", totalReads));

		
		
		
		header.add("NumCoveredPositions");
		list.add(String.format("%.0f", depthSummarizer.getNumVals()));
		header.add("CoveredPositionMeanReads");
		list.add(String.format("%.1f", depthSummarizer.getValMean()));
		header.add("GCfrac");
		list.add(String.format("%.2f", (double)this.numGCs / (double)depthSummarizer.getNumVals()));
		
		
		for (String key : result.keySet())
		{
			if (key != null)
			{
				CpgMethLevelSummarizer summarizer = result.get(key);
				header.add(String.format("pct_%s", key));
				list.add(String.format("%.3f", summarizer.getNumVals() / (double)depthSummarizer.getNumVals()));
				header.add(String.format("meth_%s", key));
				list.add(String.format("%.3f", summarizer.getValMean()));

			}
		}
		
		if (includeHeader) out.println(ListUtils.excelLine(header));
		out.println(ListUtils.excelLine(list));
		
		
	}
}
