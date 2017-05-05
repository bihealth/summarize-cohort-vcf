package de.bihealth.cubi.summarize_cohort_vcf;

import com.google.common.collect.ImmutableSet;
import de.charite.compbio.jannovar.Jannovar;
import de.charite.compbio.jannovar.pedigree.PedFileContents;
import de.charite.compbio.jannovar.pedigree.PedFileReader;
import de.charite.compbio.jannovar.pedigree.PedParseException;
import de.charite.compbio.jannovar.pedigree.PedPerson;
import de.charite.compbio.jannovar.pedigree.Sex;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

public class App {
	private static final String FOUNDER_PREFIX = "FOUNDER_";
	private static final String COHORT_AN = "COHORT_AN";
	private static final String COHORT_AC = "COHORT_AC";
	private static final String COHORT_AF = "COHORT_AF";
	private static final String COHORT_HEMI = "COHORT_Hemi";
	private static final String COHORT_HET = "COHORT_Het";
	private static final String COHORT_HOM = "COHORT_Hom";

	// Somewhat hacky but working inclusion of special chromosomes
	private static final ImmutableSet<String> X_NAMES =
			ImmutableSet.of("x", "X", "23", "chrx", "chrX", "chr23");
	private static final ImmutableSet<String> Y_NAMES =
			ImmutableSet.of("y", "Y", "24", "chry", "chrY", "chr24");
	private static final ImmutableSet<String> MT_NAMES =
			ImmutableSet.of("m", "M", "mt", "MT", "chrm", "chrM", "chrmt", "chrMT");

	ProgramOptions options = new ProgramOptions();
	PedFileContents pedFileContents;
	Set<String> allSampleNames;
	Set<String> founderSampleNames;
	Set<String> maleSampleNames = new HashSet<>();

	public void run(String[] args) {
		parseArgs(args);
		loadPedigree();
		processVcf();

		System.err.println("All done. Have a nice day!");
	}

	private void parseArgs(String[] args) {
		ArgumentParser parser = ArgumentParsers.newArgumentParser("summarize-cohort-vcf")
				.defaultHelp(true).description("Generate summarize VCF from whole-cohort VCF");
		parser.version(getVersion());
		parser.addArgument("--version").help("Show Jannovar version").action(Arguments.version());
		parser.addArgument("--input-vcf").required(true).help("Path to input VCF file");
		parser.addArgument("--output-vcf").required(true).help("Path to output VCF file");
		parser.addArgument("--pedigree").help("Optional path to pedigree file");

		Namespace ns;
		try {
			ns = parser.parseArgs(args);

			options.setPathInputVcf(ns.getString("input_vcf"));
			options.setPathInputPed(ns.getString("pedigree"));
			options.setPathOutputVcf(ns.getString("output_vcf"));
		} catch (ArgumentParserException e) {
			parser.handleError(e);
			System.exit(1);
		}
	}

	private void loadPedigree() {
		// Read pedigree and extract founder sample names if path to file given on command line
		if (options.getPathInputPed() != null) {
			final PedFileReader reader = new PedFileReader(new File(options.getPathInputPed()));
			try {
				pedFileContents = reader.read();
			} catch (PedParseException | IOException e) {
				System.err.println("Could not read pedigree. See below for technical details.");
				e.printStackTrace();
				System.exit(1);
			}

			founderSampleNames = new HashSet<>();
			allSampleNames = new HashSet<>();
			for (PedPerson pedPerson : pedFileContents.getIndividuals()) {
				allSampleNames.add(pedPerson.getName());
				if (pedPerson.isFounder()) {
					founderSampleNames.add(pedPerson.getName());
					if (pedPerson.getSex() == Sex.MALE) {
						maleSampleNames.add(pedPerson.getName());
					}
				}
			}
		}

		// Check that pedigree is consistent if read from PED file
		try (VCFFileReader reader = new VCFFileReader(new File(options.getPathInputVcf()), false)) {
			final VCFHeader fileHeader = reader.getFileHeader();

			if (allSampleNames != null) {
				Set<String> missingSampleNames = new HashSet<>(allSampleNames);
				missingSampleNames.removeAll(fileHeader.getGenotypeSamples());
				if (!missingSampleNames.isEmpty()) {
					System.err.println(
							"The following samples in the PED file but missing in the VCF files!");
					System.err.println(missingSampleNames.toString());
					System.exit(1);
				}
			} else {
				System.err.println(
						"No pedigree given, not computing founder statistics, assuming all "
								+ "samples being female and using all samples");
				allSampleNames = new HashSet<>(fileHeader.getGenotypeSamples());
			}
		}
	}

	private void processVcf() {
		System.err.println("Processing VCF file...");
		try (VCFFileReader reader = new VCFFileReader(new File(options.getPathInputVcf()), false);
				VariantContextWriter writer = constructWriter()) {
			writer.writeHeader(extendHeader(reader));

			String prevContig = null;
			for (VariantContext vc : reader) {
				if (vc.getContig().equals(prevContig)) {
					System.err.println("Starting on " + vc.getContig());
					prevContig = vc.getContig();
				}
				
				vc = processVariantContext(vc, allSampleNames, "");
				if (founderSampleNames != null) {
					vc = processVariantContext(vc, founderSampleNames, FOUNDER_PREFIX);
				}
				writer.add(vc);
			}
		}
	}

	private VariantContext processVariantContext(VariantContext vc, Set<String> sampleNames,
			String prefix) {
		if (MT_NAMES.contains(vc.getContig())) {
			return vc; // ignore MT
		}

		// Shortcuts to being X or Y chromosome
		final boolean isX = X_NAMES.contains(vc.getContig());
		final boolean isY = Y_NAMES.contains(vc.getContig());

		// Prepare variant-wise counters
		int totalChromCount = 0;
		final List<Integer> altAlleleCounts = new ArrayList<>();
		final List<Integer> altHemiCounts = new ArrayList<>();
		final List<Integer> altHetCounts = new ArrayList<>();
		final List<Integer> altHomCounts = new ArrayList<>();
		for (int i = 0; i + 1 < vc.getNAlleles(); ++i) {
			altAlleleCounts.add(0);
			altHemiCounts.add(0);
			altHetCounts.add(0);
			altHomCounts.add(0);
		}

		// Perform counting
		for (String sampleName : sampleNames) {
			final Genotype gt = vc.getGenotype(sampleName);
			final boolean isMale = maleSampleNames.contains(gt.getSampleName());
			final int chromCount = getChromCount(isMale, isX, isY);

			totalChromCount += chromCount;

			// Update counters
			//
			// Note that "1/2" will be called as heterozygous both for alternative alleles 1 and 2.
			boolean isFirst = true;
			for (Allele gtAl : gt.getAlleles()) {
				final int idx = vc.getAlleleIndex(gtAl) - 1;
				if (idx >= 0) { // is non-ref
					// increment counter for allele
					altAlleleCounts.set(idx, altAlleleCounts.get(idx) + 1);
					// increment counters for hemi, het, and hom
					if (chromCount == 1) {
						altHemiCounts.set(idx, altHemiCounts.get(idx) + 1);
						break; // don't count "1/1" twice
					} else {
						if (gt.isHet()) {
							altHetCounts.set(idx, altHetCounts.get(idx) + 1);
						} else if (gt.isHomVar()) {
							if (isFirst) { // count "1/1" as two alleles but not two hom calls
								altHomCounts.set(idx, altHomCounts.get(idx) + 1);
							}
						}
					}
				}
				isFirst = false;
			}
		}

		// Compute allele frequencies
		final List<Double> altAlleleFrequencies = new ArrayList<>();
		for (Integer count : altAlleleCounts) {
			altAlleleFrequencies.add(count.doubleValue() / totalChromCount);
		}

		// Write everything to the INFO column
		vc.getCommonInfo().putAttribute(prefix + COHORT_AN, totalChromCount);
		vc.getCommonInfo().putAttribute(prefix + COHORT_AC, altAlleleCounts);
		vc.getCommonInfo().putAttribute(prefix + COHORT_AF, altAlleleFrequencies);
		vc.getCommonInfo().putAttribute(prefix + COHORT_HEMI, altHemiCounts);
		vc.getCommonInfo().putAttribute(prefix + COHORT_HET, altHetCounts);
		vc.getCommonInfo().putAttribute(prefix + COHORT_HOM, altHomCounts);

		return vc;
	}

	private int getChromCount(boolean isMale, boolean isX, boolean isY) {
		if (isMale) {
			if (isX || isY) {
				return 1;
			} else {
				return 2;
			}
		} else {
			if (isX || !isY) {
				return 2;
			} else if (isY) {
				return 0;
			}
		}
		return 0;
	}

	private VariantContextWriter constructWriter() {
		final File outputFile = new File(options.getPathOutputVcf());
		final VariantContextWriterBuilder writerBuilder = new VariantContextWriterBuilder()
				.setCreateMD5().setOutputFile(outputFile)
				.setOptions(EnumSet.of(Options.DO_NOT_WRITE_GENOTYPES, Options.INDEX_ON_THE_FLY));
		return writerBuilder.build();
	}

	private VCFHeader extendHeader(VCFFileReader reader) {
		final VCFHeader header = reader.getFileHeader();

		String warning = "";
		if (pedFileContents == null) {
			warning = " (WARNING: no pedigree given, inaccurate results on sex chromosomes)";
		}

		String[] founderNote = { "", " (only considering founders from pedigree)" };
		String[] founderPrefix = { "", FOUNDER_PREFIX };

		for (int i = 0; i < 2; ++i) {
			final String note = founderNote[i];
			final String prefix = founderPrefix[i];

			header.addMetaDataLine(new VCFInfoHeaderLine(prefix + COHORT_AN, 1,
					VCFHeaderLineType.Integer,
					"Total number of alleles in cohort's called genotypes" + note + warning));
			header.addMetaDataLine(new VCFInfoHeaderLine(prefix + COHORT_AC, VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer, "Allele count in cohort's called genotypes" + note
							+ ", for each ALT allele, in the same order as listed" + warning));
			header.addMetaDataLine(new VCFInfoHeaderLine(prefix + COHORT_AF, VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer, "Allele Frequency in cohort, for each ALT allele"
							+ note + ", in the same order as listed" + warning));
			header.addMetaDataLine(new VCFInfoHeaderLine(prefix + COHORT_HEMI, VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer, "Cohort's Hemizygous counts, for each ALT allele"
							+ note + ", in the same order as listed" + warning));
			header.addMetaDataLine(new VCFInfoHeaderLine(prefix + COHORT_HET, VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer, "Cohort's Heterozygous counts, for each ALT allel"
							+ note + "e, in the same order as listed" + warning));
			header.addMetaDataLine(new VCFInfoHeaderLine(prefix + COHORT_HOM, VCFHeaderLineCount.A,
					VCFHeaderLineType.Integer, "Cohort's Homozygous counts, for each ALT allele"
							+ note + ", in the same order as listed" + warning));
		}

		return header;
	}

	public static void main(String[] args) {
		new App().run(args);
	}

	public static String getVersion() {
		return App.class.getPackage().getSpecificationVersion();
	}
}
