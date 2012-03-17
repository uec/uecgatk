package edu.usc.epigenome.uecgatk.pileup;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import edu.usc.epigenome.uecgatk.FractionNonidentical;

/**
 * @author benb
 *
 * This is the same idea that GATK has envisioned for Pileup2.  I am copying in their text below.
 * 
 * We also contain code to support aribitrary supergroups of read groups, which are simply Mappings
 * from a supergroup string to a list of readGroup strings.  But this makes it possible to get
 * summary information for particular read groups or read supergroups
 * 
 * We also support getting summaries by context.  How should we deal with contexts?  Do we need ot
 * get from the engine the number of upstream and downstream context positions supported?
 * 
 * 
 * From GATK:Pileup2 - notes
  - intrinsic support for samples.  Pileups are tree based data structures whose leaves are actual
pile ups of a single "sample" and whose nodes are collections of multiple sub pileups.  This will
make join and split operations very cheap.

- should be light-weight to create, and hold only minimal cached data to avoid unnecessary overhead.
Things like the number of deletions, insertions, etc shouldn't be required information.  Size will
continue to be a key cached value.  Could create a simple caching data structure that calculations lots of metrics about the pileup and was somehow
cached internally, via a "CachedRBP" structure.  This will make it very cheap and easy to filter
pileups on the fly, costing O(N) to create the filtered context.

- immutable

- support for holding neighboring reads to the left and right of the pileup itself

- unified picture for "regular" and "extended" events.  ExtendedEvents are really a special
call from the engine and have nothing to do with the data itself.

- Where should algorithms operating on the pileups go?  Two options are in the interface itself,
making it very heavy-weight but easy to access, vs. in an associated PileupOps static methods, a
la Collections.

- The Pileup2 should support in the fly filtering, so that read filters can be added at the top level
and applied at all levels of the tree.  Basically a filtering pileup would just create a new
mirrored tree with filtering applied to each node.   Very low overhead.

- Sizes could be cached as needed, so that only one pass is ever needed over the size of any pileup

- Fundamentally pileups are just collections of read-backed data.  The PileupElements contain
all of the smarts -- regular, indel, fragment-based.  We need to be able to create pileups containing
multiple subtype elements, which by necessity will need to declare their own static consensusType.  How is it
best to do this in Java?  Have a single global ENUM that enumerates all of the possible types at
compile time?  Perhaps something more dynamic?
 *
 */

// public abstract class BisulfitePileup extends Pileup2 {

public abstract class BisulfitePileup {

	
	// private variables
	Map<String,Set<String>> supergroups = null; 
	
	
	// Supergroups

	public void clearSupergroups()
	{
		supergroups = new HashMap<String,Set<String>>();
	}
	
	public void addSupergroupEntry(String supergroup, String readGroup)
	{
	}
	
	public void removeSupergroupEntry(String supergroup, String readGroup)
	{
	}
	
	
	// Methylation summaries

	
	/**
	 * @param group
	 * @return
	 */
	public FractionNonidentical getMethylation(String group)
	{
		return getMethylation(group, null);
	}
	
	/**
	 * @param group can refer to either supergroup or read group?
	 * @param context If null, return any cytosine context.
	 * @return
	 */
	
	public FractionNonidentical getMethylation(String group, String context)
	{
		FractionNonidentical out = new FractionNonidentical(0,0);
		
		return out;
	}

}
