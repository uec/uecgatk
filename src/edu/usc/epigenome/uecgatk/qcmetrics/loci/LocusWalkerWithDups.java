package edu.usc.epigenome.uecgatk.qcmetrics.loci;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.filters.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

/**
 since gatk hardcoded dup filter in, had to roll my own and remove it.
 */
@By(DataSource.READS)
@Requires({DataSource.READS,DataSource.REFERENCE, DataSource.REFERENCE_BASES})
@PartitionBy(PartitionType.INTERVAL)
@ReadFilters({UnmappedReadFilter.class,NotPrimaryAlignmentReadFilter.class,FailsVendorQualityCheckReadFilter.class})
public abstract class LocusWalkerWithDups<MapType, ReduceType> extends Walker<MapType, ReduceType> {
    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return true;    // We are keeping all the reads
    }

    // Map over the org.broadinstitute.sting.gatk.contexts.AlignmentContext
    public abstract MapType map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context);
}
