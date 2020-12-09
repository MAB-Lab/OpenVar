import os
import subprocess

class VCF:
	def __init__(self, data_dir, file_name):
		self.data_dir  = data_dir
		self.file_name = file_name
		self.file_path = os.path.join(data_dir, file_name)

	def sanitize(self, ):
		# check genome version - convert?
		# check ref-alt allele order GNOMAD
		# check 'chr' in chrome names

		# return proper feedback in case of non-compliance

	def split_per_chrom(self, ):
		# may not be necessary

	def store(self, ):
		# store sanitized vcf in db 

class OpenVar:
	def __init__(self, snpeff_path, vcf, genome_build):
		self.snpeff_jar   = os.path.join(snpeff_path, 'snpEff.jar')
		self.snpsift_jar  = os.path.join(snpeff_path, 'SnpSift.jar')
		self.snpeff_build = 'GRCh38.95_refAlt'
		self.vcf = vcf

	def run_snpeff_pipe(self, ):
		snpeff_cmd  = self.get_snpeff_cmd()
		snpeff_subproc = subprocess.Popen(snpeff_cmd.split(), shell=False, stdout=open(, 'w'))
		snpeff_subproc.wait()

		cat_cmd  = self.get_cat_cmd()
		perl_cmd = self.get_perl_cmd()
		cat_subproc  = subprocess.Popen(cat_cmd.split(), shell=False, stdout=subprocess.PIPE)
		perl_subproc = subprocess.Popen(perl_cmd.split(), shell=False, stdin=cat_subproc.stdout, stdout=subprocess.PIPE)
		cat_subproc.stdout.close()

		snpsift_cmd  = self.get_snpsift_cmd()
		snpsift_subproc = subprocess.Popen(snpsift_cmd, shell=True, stdin=perl_subproc.stdout, stdout=open(chr_cmds[chrom]['annOnePerLine'], "w"))

	    perl_subproc.stdout.close()
	    snpsift_subproc.wait()		

	def get_snpeff_cmd(self, ):
		cmd = 'java -Xmx12g -jar {snpeff_jar} -v {build} {vcf_path}'.format(
				snpeff_jar = self.snpeff_jar,
				build = self.snpeff_build,
				vcf_path = self.vcf.file_path
			)
		return cmd

	def get_cat_cmd(self, ):
		cmd = 'cat {output_vcf}'.format(output_vcf=self.vcf.file_name.replace('.vcf', '.ann.vcf'))
		return cmd

	def get_perl_cmd(self, ):
		cmd = os.path.join(self.snpeff_path, 'scripts', 'vcfEffOnePerLine.pl')
		return cmd

	def get_snpsift_cmd(self, ):
		snpsift_fields = 'CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].HGVS_P" "ANN[*].FEATUREID" "ANN[*].GENE" "ANN[*].ERRORS" "ANN[*].HGVS_C"'
		cmd ='java -jar {snpsift_jar} extractFields - {snpsift_fields}'.format(
				snpsift_jar = self.snpsift_jar,
				snpsift_fields = snpsift_fields
			)
		return cmd



class OPVReport:
	def __init__(self, opv):

	def write_tabular(self, ):

	def write_html(self, ):