package edu.usc.epigenome.uecgatk.benWalkers;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.bio.seq.StrandedFeature.Strand;

import edu.usc.epigenome.uecgatk.FractionNonidentical;
import edu.usc.epigenome.uecgatk.IupacPatterns;

public class ReadWithCpgMethsQuadrants extends ReadWithCpgMeths {

	protected static final FractionNonidentical all = new FractionNonidentical(100,100); 
	protected static final FractionNonidentical none = new FractionNonidentical(0,0); 
	
    private static final Map<String, double[]> contextEdgesMap = new HashMap<String, double[]>();
    static {
    	contextEdgesMap.put("WCG", new double[] { 0.0, 0.1, 0.2, 0.5, 0.8, 1.0 });
    	contextEdgesMap.put("HCG", new double[] { 0.0, 0.1, 0.2, 0.5, 0.8, 1.0 });
       	contextEdgesMap.put("CCG", new double[] { 0.0, 0.1, 0.2, 0.5, 0.8, 1.0 });
       	contextEdgesMap.put("GCH", new double[] { 0.0, 0.1, 0.9, 1.0 });
    }
	
	public ReadWithCpgMethsQuadrants(Strand inStrand, String inChrom) {
		super(inStrand, inChrom);
	}

	
	public static Map<String, FractionNonidentical> methLevelsFractions(ReadWithCpgMeths inRead, IupacPatterns patterns) 
	{
//		return methLevelsFractions(inRead, patterns, 0.1, 0.9);
		return methLevelsFractions(inRead, patterns, contextEdgesMap);
	}

	public static Map<String, FractionNonidentical> methLevelsFractions(ReadWithCpgMeths inRead, IupacPatterns patterns,double lowMax, double highMin)
	{
		Map<String, FractionNonidentical> out = new HashMap<String, FractionNonidentical>();
		Map<String, FractionNonidentical> origLevels =  inRead.methLevelsFractions(patterns);
		out.putAll(origLevels);
		List<String> contexts = new ArrayList(origLevels.keySet());

		for (int i = 0; i < (contexts.size()-1); i++)
		{
			String con_i = contexts.get(i);
			double level_i = origLevels.get(con_i).doubleValue();
			
			for (int j = (i+1); j < contexts.size(); j++)
			{
				String con_j = contexts.get(j);
				double level_j = origLevels.get(con_j).doubleValue();
				
				String con_new = null;
				final String formatStr = "%s%s%.2f-%s%s%.2f";

				con_new = String.format(formatStr,con_i,"lt",lowMax,con_j,"lt",lowMax);
				out.put(con_new, ((level_i < lowMax) && (level_j < lowMax)) ? all : none);

				con_new = String.format(formatStr,con_i,"lt",lowMax,con_j,"gt",highMin);
				out.put(con_new, ((level_i < lowMax) && (level_j > highMin)) ? all : none);

				con_new = String.format(formatStr,con_i,"gt",highMin,con_j,"lt",lowMax);
				out.put(con_new, ((level_i > highMin) && (level_j < lowMax)) ? all : none);
				
				con_new = String.format(formatStr,con_i,"gt",highMin,con_j,"gt",highMin);
				out.put(con_new, ((level_i > highMin) && (level_j > highMin)) ? all : none);
			}
		}
		return out;
	}

	/**
	 * 
	 * @param inRead
	 * @param patterns
	 * @param edgesByContext  edges for each "context" in context map.
	 * @return
	 */
	public static Map<String, FractionNonidentical> methLevelsFractions(ReadWithCpgMeths inRead, IupacPatterns patterns, Map<String, double[]> edgesByContext)
	{
		Map<String, FractionNonidentical> out = new HashMap<String, FractionNonidentical>();
		Map<String, FractionNonidentical> origLevels =  inRead.methLevelsFractions(patterns);
		out.putAll(origLevels);
		List<String> contexts = new ArrayList(origLevels.keySet());

		for (int i = 0; i < (contexts.size()-1); i++)
		{
			String con_i = contexts.get(i);
			double level_i = origLevels.get(con_i).doubleValue();
			double[] edges_i = edgesByContext.get(con_i);

			for (int j = (i+1); j < contexts.size(); j++)
			{
				String con_j = contexts.get(j);
				double level_j = origLevels.get(con_j).doubleValue();
				double[] edges_j = edgesByContext.get(con_j);

				String con_new = null;
				final String formatStr = "%sgt%.2flt%.2f-%sgt%.2flt%.2f";
				
				
				for (int edgenumi = 0; edgenumi < edges_i.length-1; edgenumi++)
				{
					for (int edgenumj = 0; edgenumj < edges_j.length-1; edgenumj++)
					{
						double ei1 = edges_i[edgenumi];
						double ei2 = edges_i[edgenumi+1];
						double ej1 = edges_j[edgenumj];
						double ej2 = edges_j[edgenumj+1];

						con_new = String.format(formatStr,con_i,ei1,ei2,con_j,ej1,ej2);
						boolean passes = true;
						// First context
						passes &= (level_i <= ei2);
						if (edgenumi == 0) // first
						{
							passes &= (level_i >= ei1);
						}
						else
						{
							passes &= (level_i > ei1);
						}
						// Second context
						passes &= (level_j <= ej2);
						if (edgenumj == 0) // first
						{
							passes &= (level_j >= ej1);
						}
						else
						{
							passes &= (level_j > ej1);
						}
						
						out.put(con_new, (passes) ? all : none); 
					}
				}

			}
		}

		return out;
	}

	
	
}
