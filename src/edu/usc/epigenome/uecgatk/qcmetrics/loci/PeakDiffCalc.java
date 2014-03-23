package edu.usc.epigenome.uecgatk.qcmetrics.loci;

import java.util.List;

public interface PeakDiffCalc
{
	public float calcDiff(List<SinglePeak> peaks);
}
