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
	

	
	public ReadWithCpgMethsQuadrants(Strand inStrand, String inChrom) {
		super(inStrand, inChrom);
	}

	
	public static Map<String, FractionNonidentical> methLevelsFractions(ReadWithCpgMeths inRead, IupacPatterns patterns) 
	{
		return methLevelsFractions(inRead, patterns, 0.1, 0.9);
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

	
	
}
