package edu.usc.epigenome.uecgatk;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broadinstitute.sting.utils.BaseUtils;

import edu.usc.epigenome.genomeLibs.ListUtils;

public class IupacPatterns {

//	protected List<String> iupacs = new ArrayList<String>();

	protected Map<String, String> map = new HashMap<String,String>(); // seq -> iupac
	

//	/**
//	 * @param iupac
//	 * If an example matches multiple, we return the first one registered.
//	 */
//	public void register(String iupac)
//	{
//		iupacs.add(iupac);
//	}
//		
//	/**
//	 * @return the iupacs
//	 */
//	public List<String> getIupacs() {
//		return iupacs;
//	}

	/**
	 * @param iupac
	 * If an example matches multiple, we return the first one registered.
	 */
	public void register(String iupac)
	{
		// Add all permutations to map.
		List<String> perms = getPermutations("",iupac);
		System.err.printf("Found %d permutations for %s:  %s\n", perms.size(), iupac, ListUtils.excelLine(perms));

		for (String perm : perms)
		{
			// Don't add if a permuation is already registered to another iupac
			if (!this.map.containsKey(perm))
			{
				this.map.put(perm, iupac);
			}
		}
	}

	public Collection<String> getRegistered()
	{
		return this.map.values();
	}

	public String firstMatch(String seq)
	{
		String out = map.get(seq.toUpperCase());
		return out;
	}



	
	
	
	
	protected static List<String> getPermutations(String prefix, String iupac)
	{
		List<String> out = new ArrayList<String>();

		// End of recursion
		if (iupac.length()==0)
		{
			out.add(prefix);
		}
		else
		{
			char iupacChar = iupac.charAt(0);
			String iupacRest = (iupac.length()==1) ? "" : iupac.substring(1);

			// Get bare suffixes
			List<String> suffixes = getPermutations("",iupacRest);

			for (char c : getIupacBases(iupacChar))
			{
				for (String suf : suffixes)
				{
					out.add(String.format("%s%c%s", prefix, c, suf));
				}
			}
		}
		
		return out;
	}
	
	
	public static char[] getIupacBases(char in)
	{
		char[] out = new char[]{};
		
		switch (Character.toUpperCase(in))
		{
		case 'A':
		case 'C':
		case 'T':
		case 'G':
			out = new char[]{in};
			break;
		case 'N':
			out = new char[]{'A','C','T','G'};
			break;
		case 'Y':
			out = new char[]{'C','T'};
			break;
		case 'U':
			out = new char[]{'A','G'};
			break;
		case 'W':
			out = new char[]{'A','T'};
			break;
		case 'S':
			out = new char[]{'C','G'};
			break;
		case 'H':
			out = new char[]{'C','T','A'};
			break;
		}
		
		return out;
	}
	
	
}
