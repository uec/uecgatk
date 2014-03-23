/*
 * Copyright (c) 2012 USC Epigenome Center
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package edu.usc.epigenome.uecgatk.filters;

import net.sf.samtools.SAMRecord;

import org.broadinstitute.sting.gatk.filters.ReadFilter;

/**
 * Filter out reads that are improperly paired.
 *
 * @author zack
 * @version 0.1
 */

public class NumberHitsEqualsOneFilter extends ReadFilter
{

	@Override
	public boolean filterOut(SAMRecord read)
	{
		try
		{
			return (read.getNotPrimaryAlignmentFlag() || read.getIntegerAttribute("NH") != 1);
		}
		catch (Exception e)
		{
			return false;
		}
	}
}