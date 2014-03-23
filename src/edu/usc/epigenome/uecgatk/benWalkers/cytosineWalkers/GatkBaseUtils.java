package edu.usc.epigenome.uecgatk.benWalkers.cytosineWalkers;

import org.broadinstitute.sting.utils.BaseUtils;

public class GatkBaseUtils {
	public static char CharFromBaseByte(byte inBase)
	{
		return (char) BaseUtils.baseIndexToSimpleBase(BaseUtils.extendedBaseToBaseIndex(inBase));
	}
	
	public static byte BaseByteFromChar(char inBaseChar)
	{
		byte out = BaseUtils.Base.N.base;
		switch (inBaseChar)
		{
		case 'A':
		case 'a':
			out = BaseUtils.Base.A.base;
			break;
		case 'C':
		case 'c':
			out = BaseUtils.Base.C.base;
			break;
		case 'T':
		case 't':
			out = BaseUtils.Base.T.base;
			break;
		case 'G':
		case 'g':
			out = BaseUtils.Base.G.base;
			break;
		case 'N':
		case 'n':
		case '0':
		default:
			out = BaseUtils.Base.N.base;
			break;
		}
		
		return out;
	}
	
    /**
     * Copied from org.broadinstitute.sting.utils.BaseUtils::mostFrequentBaseFraction
     * Finds the most frequent base in the sequence
     *
     * @param sequence  the read sequence
     * @return  the percentage of the read that's made up of the most frequent base
     */
    static public byte mostFrequentBase(byte[] sequence) {
        int[] baseCounts = new int[4];

        for ( byte base : sequence ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);

            if (baseIndex >= 0) {
                baseCounts[baseIndex]++;
            }
        }

        int mostFrequentBaseIndex = 0;
        for (int baseIndex = 1; baseIndex < 4; baseIndex++) {
            if (baseCounts[baseIndex] > baseCounts[mostFrequentBaseIndex]) {
                mostFrequentBaseIndex = baseIndex;
            }
        }

        return BaseUtils.baseIndexToSimpleBase(mostFrequentBaseIndex);
    }
}
