package BuildConsensus;

import java.io.*;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import jaligner.Sequence;
import jaligner.matrix.MatrixLoaderException;

public class BuildMain {

	private static File bamFile;
	private static File bamIndex;
	private static File sequenceOut;
	private static File referenceSequence;
	private static File referenceIndex;
	private static String chr;
	private static int start;
	private static int end;
	
	public static void main(String[] args) throws IOException, InterruptedException, MatrixLoaderException, URISyntaxException {

		loadArgs(args);
		
		IndexedFastaSequenceFile reference = new IndexedFastaSequenceFile(referenceSequence, new FastaSequenceIndex(referenceIndex));
		
		BufferedWriter output = new BufferedWriter(new FileWriter(sequenceOut));
		
//		Sequence meiLoaded = loadSequence(meiSequence);
		Sequence seqLoaded = loadSequence(chr, start, end, reference);
		List<SAMRecord> records = loadRecords(bamFile, bamIndex, chr, start, end);
		String sequence = BuildSpecies.getSpecies(seqLoaded, records, start, end);
		output.write(">chr" + chr + ":" + start + "-" + end + "\n" + sequence);
		output.newLine();
		output.flush();
		output.close();
		
	}

	private static Sequence loadSequence (String chr, int start, int end, IndexedFastaSequenceFile reference) {
		
		return new Sequence(reference.getSubsequenceAt(chr, start, end).getBaseString());
		
	}
	private static List<SAMRecord> loadRecords (File bam, File index, String chr, int start, int stop) throws IOException {
		
		SamReader bamReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(bam).index(index));
		SAMRecordIterator itr = bamReader.queryOverlapping(chr, start, stop);
		List<SAMRecord> records = new ArrayList<SAMRecord>();
		
		while (itr.hasNext()) {
			
			SAMRecord next = itr.next();
			if (next.getReadUnmappedFlag() == false) {
				
				records.add(next);
				
			}
		
		}
		
		bamReader.close();
		return records;
		
	}
	
	private static void loadArgs (String args[]) {
		
		Options options = new Options();
		List<Option> reqOptions = new ArrayList<Option>();
		reqOptions.add(new Option("bamfile", true, "Bam file to get consensus seq from (needs .bai index)."));
		reqOptions.add(new Option("o", true, "FASTA format output of consensus seq."));
		reqOptions.add(new Option("r", true, "Reference file bam is aligned to (needs .fai index)."));
		reqOptions.add(new Option("chr", true, "Chromosome."));
		reqOptions.add(new Option("start", true, "Start coordinate."));
		reqOptions.add(new Option("end", true, "End coordinate."));
		
		for (Option opt : reqOptions) {
			opt.setRequired(true);
			options.addOption(opt);
		}
		
		options.addOption(new Option("maj", true, "Percent majority to call consensus [0.50]"));
		
		CommandLineParser parser = new BasicParser();
		CommandLine cmd = null;
		
		try {
			 cmd = parser.parse(options, args);
		} catch (org.apache.commons.cli.ParseException e) {
			System.out.println();
			ThrowHelp(e.getMessage(), options);
		}
		
		bamFile = new File(cmd.getOptionValue("bamfile"));
		bamIndex = new File(cmd.getOptionValue("bamfile") + ".bai");
		sequenceOut = new File(cmd.getOptionValue("o"));
		referenceSequence = new File(cmd.getOptionValue("r"));
		referenceIndex = new File(cmd.getOptionValue("r") + ".fai");
		chr = cmd.getOptionValue("chr");
		start = Integer.parseInt(cmd.getOptionValue("start"));
		end = Integer.parseInt(cmd.getOptionValue("end"));
		
	}
	
	private static void ThrowHelp(String message, Options options) {
		System.out.println();
		System.out.println(message);
		HelpFormatter formatter = new HelpFormatter();
		System.out.println();
		formatter.printHelp(9999, "java -jar BuildConsensus.jar <options>", "BuildConsensus - BuildConsensus sequences using majority rule from aligned bam files", options, "\n\n-help will print this message and exit");
		System.exit(1);
	}
	
}
