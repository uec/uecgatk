package edu.usc.epigenome.uecgatk.benMiscScripts;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

import org.broadinstitute.sting.utils.collections.Pair;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import BisulfiteCytosines.CpgPair;

import edu.usc.epigenome.genomeLibs.MatUtils;

public class KaplanNucleosomePullStrongestRegions {

	private static String C_USAGE = "Use: KaplanNucleosomePullStrongestRegions -minQuantile 85 file1.txt ";

	@Option(name="-minQuantile",usage="1-100 (default 15)")
	private int minQuantile = 85;
	@Option(name="-minScore",usage="0.0-1.0 (default false, calculate from quantile)")
	private double minScore = Double.NaN;
	@Option(name="-dontWriteFile",usage="(default false)")
	private boolean dontWriteFile = false;
	@Option(name="-nucWindSize",usage="don't output positions within this distance of each other (default=185)")
	private int nucWindSize = 185;
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

	protected void outputRegions(String fn, double inMinScore) 
	throws Exception
	{
		long[] scores = processFile(fn, inMinScore, true);
	}

	protected double determineMinScore(String fn, int inMinQuantile) 
	throws Exception
	{
		long[] scores = processFile(fn, 0.0, false);
		
		// Work backwards until we're below the quantile
		long total = 0;
		int i;
		for (i=0; i < scores.length; i++) total+=scores[i];
		long minNumNeeded = (long)Math.round((double)total*(double)inMinQuantile/100.0);
		System.err.printf("Find %d to include in %dth percentile (%d total)\n", minNumNeeded, inMinQuantile, total);
		long totalIncluded = 0;
		for (i = scores.length-1; (totalIncluded<minNumNeeded) && (i>=0); i--)
		{
			totalIncluded += scores[i];
			//System.err.printf("\ti=%d, totalIncluded=%d\n", i,totalIncluded);
		}
		
		i++;
		double minScore = (double)i / (double)(scores.length);
		return minScore;
	}
	
	protected long[] processFile(String fn, double inMinScore, boolean writeFile)
	throws Exception
	{
		// Now read the file
				
		boolean gzipped = fn.endsWith(".gz");
		BufferedReader reader = (gzipped) ? new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(fn))))
				: new BufferedReader(new FileReader(fn));

		long[] out = new long[1000];
		
		String line = null;
		int lineCount = 0;
		WindowTracker windTrack = new WindowTracker(this.nucWindSize);
		double lastScoreOutput = 0.0;
		int lastScoreOutputCoord = 1;
		String lastChr = "";
		while ((line = reader.readLine()) != null)
		{
			lineCount++;
			if ((lineCount%1000)==0) System.err.printf("On line %d...\n",lineCount);
			
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
			
			if (!chr.equals(lastChr))
			{
				// New chrom
				windTrack.init();
			}
			lastChr = chr;

			for (int coord = start; coord <= end; coord++)
			{
				//				// *****Debugging
				//				if (writeFile && (coord > 2000))
				//				{
				//					reader.close();
				//					return out;
				//				}
				//				// *****Debugging
						
				
				
				if (!tok.hasMoreTokens())
				{
					System.err.printf("Why does line %d not have a tokent for coord %d?\n", lineCount, coord);
				}
				else
				{
					String token = tok.nextToken(";");
					double score = Double.parseDouble(token);

					// Add to counts
					int index = Math.min(999,(int)(score*1000.0));
					out[index]++;

					if (writeFile)
					{
						// Adjust window info
						windTrack.addToWind(coord, score);

						//					if (score >= inMinScore)
						//if (bestScoreInWindCoord==(coord-inNucWindSize))
						
						int distToLast = (windTrack.getBestInWindowCoord()-lastScoreOutputCoord);
						boolean bestScoreFarEnough = (distToLast > this.nucWindSize);
						// System.err.printf("\tchr%s, %d, %.4f (best=%d, %.4f, dist to last=%d)\n",chr,coord,score,windTrack.getBestInWindowCoord(), windTrack.getBestInWindowScore(), distToLast);
						if (bestScoreFarEnough)
						{
							//System.err.printf("Best in window %d, %.3f (dist to last=%d)\n", windTrack.getBestInWindowCoord(), windTrack.getBestInWindowScore(), distToLast);
							boolean bestScoreAboveMin = (windTrack.getBestInWindowScore()>=inMinScore); 
							if (bestScoreAboveMin) 
							{
								//System.err.printf("\tOutputting best score\n");
								//	chr1    BSPP    exon    149505458       149505467       .       +       .
								int bestIndex = Math.min(999,(int)(windTrack.getBestInWindowScore()*1000.0));
								System.out.printf("chr%s\tKap08\texon\t%d\t%d\t%d\t+\t.\n", chr, windTrack.getBestInWindowCoord(), windTrack.getBestInWindowCoord(), bestIndex);
								
								lastScoreOutputCoord = windTrack.getBestInWindowCoord();
								lastScoreOutput = windTrack.getBestInWindowScore();
							}
						}
					}
				}
			}
		}

		reader.close();
		return out;
	}
	
	
	public class WindowTracker
	{
		protected int windLength = 0;
		
		public LinkedList<Pair<Integer,Double>> windowScores;
		int bestCoordInWind;
		double bestScoreInWind;


		public WindowTracker(int inWindLength) {
			super();
			this.windLength = inWindLength;
			init();
		}

		private void init() 
		{
			windowScores = new LinkedList<Pair<Integer,Double>>();
			bestCoordInWind = -1;
			bestScoreInWind = Double.NEGATIVE_INFINITY;
		}
		
		
		/**
		 * @return the bestInWindowCoord
		 */
		public int getBestInWindowCoord() {
			return bestCoordInWind;
		}


		/**
		 * @return the bestInWindowScore
		 */
		public double getBestInWindowScore() {
			return bestScoreInWind;
		}


		private void addToWind(int inCoord, double inScore)
		{
			// See if we have the best score
			if (inScore > bestScoreInWind)
			{
				setBestScore(inCoord,inScore);
			}
			
			// Eat off the head if necessary. Check if we are eating the best score.
			if (windowScores.size()>=this.windLength)
			{
				Pair<Integer,Double> removePair = windowScores.remove();
				if (removePair.getFirst() == bestCoordInWind)
				{
					resetBestScore();
				}
			}
			
			
			// Then add the new one
			Pair<Integer,Double> pair = new Pair<Integer,Double>(inCoord,inScore);
			windowScores.add(pair);
			
			// Debugging
			//System.err.println("\t\t" + desc());

		}
		
		private void resetBestScore()
		{
			//System.err.println("resetBestScore()\n\t" + desc());
			bestCoordInWind = -1;
			bestScoreInWind = Double.NEGATIVE_INFINITY;
			for (Pair<Integer,Double> cur : windowScores)
			{
				if (cur.getSecond() > bestScoreInWind)
				{
					setBestScore(cur.first, cur.second);
				}
			}
		}
		
		private void setBestScore(int inCoord, double inScore)
		{
			bestCoordInWind = inCoord;
			bestScoreInWind = inScore;
		}
		
		public String desc()
		{
			String out = String.format("windTracker:%d els (%d-%d) bestCoord=%d  bestScore=%.3f", 
					windowScores.size(), windowScores.getFirst().first, windowScores.getLast().first, bestCoordInWind, bestScoreInWind);
			return out;
		}
		
	}
	

}


