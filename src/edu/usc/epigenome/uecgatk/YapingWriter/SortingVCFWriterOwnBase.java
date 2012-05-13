package edu.usc.epigenome.uecgatk.YapingWriter;

import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.PriorityBlockingQueue;

import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;

import edu.usc.epigenome.uecgatk.BisSNP.BaseUtilsMore;



public abstract class SortingVCFWriterOwnBase implements VCFWriter {

	// The VCFWriter to which to actually write the sorted VCF records
    private VCFWriter innerWriter = null;

    // the current queue of un-emitted records
    private PriorityBlockingQueue<VCFRecord> queue = null;

    // The locus until which we are permitted to write out (inclusive)
    protected Integer mostUpstreamWritableLoc;
    protected String cachedCurrentChr = null;
    
    protected boolean enableDiscreteLoci = false;
    
    
    protected static final int BEFORE_MOST_UPSTREAM_LOC = 0; // No real locus index is <= 0

    // The set of chromosomes already passed over and to which it is forbidden to return
    private Set<String> finishedChromosomes = null;

    // Should we call innerWriter.close() in close()
    private boolean takeOwnershipOfInner;

    /**
     * create a local-sorting VCF writer, given an inner VCF writer to write to
     *
     * @param innerWriter        the VCFWriter to write to
     * @param takeOwnershipOfInner Should this Writer close innerWriter when it's done with it
     */
    public SortingVCFWriterOwnBase(VCFWriter innerWriter, boolean takeOwnershipOfInner) {
        this.innerWriter = innerWriter;
        this.queue = new PriorityBlockingQueue<VCFRecord>(100000000, new VariantContextComparator());
        this.mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
        this.finishedChromosomes = new TreeSet<String>();
        this.takeOwnershipOfInner = takeOwnershipOfInner;
    }

    public SortingVCFWriterOwnBase(VCFWriter innerWriter) {
        this(innerWriter, false); // by default, don't own inner
    }

    public void writeHeader(VCFHeader header) {
        innerWriter.writeHeader(header);
    }

    /**
     * attempt to close the VCF file; we need to flush the queue first
     */
    public void close() {
        stopWaitingToSort();

        if (takeOwnershipOfInner)
            innerWriter.close();
    }

    private void stopWaitingToSort() {
        emitRecords(true);
        mostUpstreamWritableLoc = BEFORE_MOST_UPSTREAM_LOC;
    }

    public void emitSafeRecords() {
        emitRecords(false);
    }

    protected void noteCurrentRecord(VariantContext vc) {
        // did the user break the contract by giving a record too late?
        if (mostUpstreamWritableLoc != null && (cachedCurrentChr != null && cachedCurrentChr.equalsIgnoreCase(vc.getChr()) && vc.getStart() < mostUpstreamWritableLoc) && !enableDiscreteLoci) // went too far back, since may have already written anything that is <= mostUpstreamWritableLoc
            throw new IllegalArgumentException("Permitted to write any record upstream of position : " + cachedCurrentChr + ": " + mostUpstreamWritableLoc + ", but a record at " + vc.getChr() + ":" + vc.getStart() + " was just added.");
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     * @param refBase the ref base
     */
    public synchronized void add(VariantContext vc) {
        /* Note that the code below does not prevent the successive add()-ing of: (chr1, 10), (chr20, 200), (chr15, 100)
           since there is no implicit ordering of chromosomes:
         */
    	if(cachedCurrentChr == null){
    		cachedCurrentChr = vc.getChr();
    	}
        VCFRecord firstRec = queue.peek();
        if (firstRec != null && !vc.getChr().equals(firstRec.vc.getChr())) { // if we hit a new contig, flush the queue
           // if (finishedChromosomes.contains(vc.getChr()) && !enableDiscreteLoci)
           //     throw new IllegalArgumentException("Added a record at " + vc.getChr() + ":" + vc.getStart() + ", but already finished with chromosome" + vc.getChr());

            finishedChromosomes.add(firstRec.vc.getChr());
            cachedCurrentChr = vc.getChr();
            stopWaitingToSort();
        }

        noteCurrentRecord(vc); // possibly overwritten

        queue.add(new VCFRecord(vc));
        emitSafeRecords();
    }
    
    

    private void emitRecords(boolean emitUnsafe) {
    	if(mostUpstreamWritableLoc != 0 && mostUpstreamWritableLoc % 100000 ==0)
    		System.err.println(mostUpstreamWritableLoc + "\t" + queue.peek().vc.getStart());
    	while (!queue.isEmpty()) {
            VCFRecord firstRec = queue.peek();

            // No need to wait, waiting for nothing, or before what we're waiting for:
            if (emitUnsafe || mostUpstreamWritableLoc == null || firstRec.vc.getStart() <= mostUpstreamWritableLoc || !cachedCurrentChr.equalsIgnoreCase(firstRec.vc.getChr()) || enableDiscreteLoci) {
                queue.poll();
                innerWriter.add(firstRec.vc);
            }
            else {
                break;
            }
        }
    }

    /**
     * Gets a string representation of this object.
     * @return
     */
    @Override
    public String toString() {
        return getClass().getName();
    }

    private static class VariantContextComparator implements Comparator<VCFRecord> {
        public int compare(VCFRecord r1, VCFRecord r2) {
            if(!r1.vc.getChr().equalsIgnoreCase(r2.vc.getChr()))
            	return BaseUtilsMore.contigToIndex(r1.vc.getChr())-BaseUtilsMore.contigToIndex(r1.vc.getChr());
        	return r1.vc.getStart() - r2.vc.getStart();
        }
    }

    private static class VCFRecord {
        public VariantContext vc;
       // public byte refBase;

        public VCFRecord(VariantContext vc) {
            this.vc = vc;
          //  this.refBase = refBase;
        }
    }
  
}
