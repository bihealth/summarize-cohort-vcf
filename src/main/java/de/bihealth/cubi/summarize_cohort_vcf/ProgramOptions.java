package de.bihealth.cubi.summarize_cohort_vcf;

/**
 * Program command line options
 * 
 * @author <a href="mailto:manuel.holtgrewe@bihealth.de">Manuel Holtgrewe</a>
 */
public class ProgramOptions {

	private String pathInputVcf;
	private String pathInputPed;
	private String pathOutputVcf;
	private boolean writeSamples = false;

	public ProgramOptions() {
	}

	public String getPathInputVcf() {
		return pathInputVcf;
	}

	public void setPathInputVcf(String pathInputVcf) {
		this.pathInputVcf = pathInputVcf;
	}

	public String getPathInputPed() {
		return pathInputPed;
	}

	public void setPathInputPed(String pathInputPed) {
		this.pathInputPed = pathInputPed;
	}

	public String getPathOutputVcf() {
		return pathOutputVcf;
	}

	public void setPathOutputVcf(String pathOutputVcf) {
		this.pathOutputVcf = pathOutputVcf;
	}

	public boolean isWriteSamples() {
		return writeSamples;
	}

	public void setWriteSamples(boolean writeSamples) {
		this.writeSamples = writeSamples;
	}

	@Override
	public String toString() {
		return "ProgramOptions [pathInputVcf=" + pathInputVcf + ", pathInputPed=" + pathInputPed
				+ ", pathOutputVcf=" + pathOutputVcf + ", writeSamples=" + writeSamples + "]";
	}

}
