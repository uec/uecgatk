package edu.usc.epigenome.uecgatk.benMiscScripts;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.usc.epigenome.genomeLibs.MatUtils;

public class KaplanNucleosomePullStrongestRegions {

	private static String C_USAGE = "Use: KaplanNucleosomePullStrongestRegions -minQuantile 85 file1.txt ";

	@Option(name="-minQuantile",usage="1-100 (default 15)")
	private int minQuantile = 85;
	@Option(name="-minScore",usage="0.0-1.0 (default false, calculate from quantile)")
	private double minScore = Double.NaN;
	@Option(name="-dontWriteFile",usage="(default false)")
	private boolean dontWriteFile = false;
	// receives other command line parameters than options
	@Argument
	private List<String> arguments = new ArrayList<String>();


	public static void main(String[] args)
	throws Exception
	{
		 new KaplanNucleosomePullStrongestRegions().doMain(args);
	}

	public void doMain(String[] args)
	throws Exception
	{
		CmdLineParser parser = new CmdLineParser(this);
		// if you have a wider console, you could increase the value;
		// here 80 is also the default
		parser.setUsageWidth(80);
		try
		{
			parser.parseArgument(args);

			if(arguments.size() != 1) {
				System.err.println(C_USAGE);
				System.exit(1);
			}

		}
		catch (CmdLineException e)
		{
			System.err.println(e.getMessage());
			System.err.println(C_USAGE);
			// print the list of available options
			parser.printUsage(System.err);
			System.err.println();
			return;
		}	



		
		String fn = arguments.get(0);
		
		if (Double.isNaN(minScore))
		{
			System.err.printf("On pass 1, getting %dth quantile cutoff\n", this.minQuantile);
			minScore = determineMinScore(fn, this.minQuantile);
		}
		System.err.printf("\tminScore=%.3f\n", minScore);
		if (!this.dontWriteFile)
		{
			System.err.printf("On pass 2, outputting gtf file\n", this.minQuantile);
			outputRegions(fn, minScore);
		}


	}

	static protected void outputRegions(String fn, double inMinScore) 
	throws Exception
	{
		int[] scores = processFile(fn, inMinScore, true);
		
		
	}

	static protected double determineMinScore(String fn, int inMinQuantile) 
	throws Exception
	{
		int[] scores = processFile(fn, 0.0, false);
		
		// Work backwards until we're below the quantile
		int total = MatUtils.nanSum(scores);
		int minNumNeeded = (int)Math.round((double)total*(double)inMinQuantile/100.0);
		System.err.printf("Find %d to include in %dth percentile (%d total)\n", minNumNeeded, inMinQuantile, total);
		int totalIncluded = 0;
		int i = 0;
		for (i = scores.length-1; (totalIncluded<minNumNeeded) && (i>=0); i--)
		{
			totalIncluded += scores[i];
			//System.err.printf("\ti=%d, totalIncluded=%d\n", i,totalIncluded);
		}
		
		i++;
		double minScore = (double)i / (double)(scores.length);
		return minScore;
	}
	
	static protected int[] processFile(String fn, double inMinScore, boolean writeFile)
	throws Exception
	{
		// Now read the file
				
		boolean gzipped = fn.endsWith(".gz");
		BufferedReader reader = (gzipped) ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fn))))
				: new BufferedReader(new FileReader(fn));

		int[] out = new int[1000];
		
		String line = null;
		int lineCount = 0;
		while ((line = reader.readLine()) != null)
		{
			lineCount++;
			if ((lineCount%100000)==0) System.err.printf("On line %d...\n",lineCount);
			
//			String[] columns = line.split("\t");
//			if (columns.length!=4)
//			{
//				System.err.printf("Why does line %d have %d columns?\n",lineCount, columns.length);
//			}
//			else
//			{
//				String chr = columns[0];
//				int start = Integer.parseInt(columns[1]);
//				int end = Integer.parseInt(columns[2]);
//				StringTokenizer tok = new StringTokenizer(columns[3],";");
			
			StringTokenizer tok = new StringTokenizer(line);
			String chr = tok.nextToken("\t");
			int start = Integer.parseInt(tok.nextToken("\t"));
			int end = Integer.parseInt(tok.nextToken("\t"));

			for (int coord = start; coord <= end; coord++)
			{
				//System.err.printf("chr%s, %d\n",chr,coord);
				if (!tok.hasMoreTokens())
				{
					System.err.printf("Why does line %d not have a tokent for coord %d?\n", lineCount, coord);
				}
				else
				{
					String token = tok.nextToken(";");
					double score = Double.parseDouble(token);
					//if (coord == 200) System.err.printf("chr%s, %d\tval=%.4f\n",chr,coord,score);
					
					int index = Math.min(999,(int)(score*1000.0));
					out[index]++;
					
					if (writeFile && (score >= inMinScore))
					{
						//	chr1    BSPP    exon    149505458       149505467       .       +       .
						System.out.printf("chr%s\tKap08\texon\t%d\t%d\t%d\t+\t.\n", chr, coord, coord, index);
					}
				}
			}
		}
		
		reader.close();
		return out;
	}

}


