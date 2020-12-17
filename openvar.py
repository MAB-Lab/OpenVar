import os
import csv
import subprocess
import multiprocessing
import itertools as itt

chrom_names = set([str(x) for x in range(1,23)] + ['X', 'Y', 'MT'])

class VCF:
	def __init__(self, data_dir, file_name, study_name):
		self.data_dir   = data_dir
		if not os.path.exists(self.data_dir):
			try:
				os.mkdir(self.data_dir)
			except:
				print("can't create data dir...")

		self.vcf_splits_dir = os.path.join(self.data_dir, 'vcf_splits')
		if not os.path.exists(self.vcf_splits_dir):
			os.mkdir(self.vcf_splits_dir)

		self.file_name  = file_name
		self.file_path  = os.path.join(data_dir, file_name)
		self.study_name = study_name
		self.parse_vcf()

	def parse_vcf(self):
		vcf_ls = []
		with open(self.file_path, 'r') as f:
			reader = csv.reader(f, delimiter='\t')
			for n, row in emuerate(reader):
				vcf_ls.append(row)
		self.vcf_ls = sorted(vcf_ls, key=lambda x: x[0])

	def check(self):
		# check genome version - convert?
		# check ref-alt allele order GNOMAD
		# check 'chr' in chrome names

		# return proper feedback in case of non-compliance
		pass

	def check_alt_ref_alleles(self):
		pass
		
	def convert_hg19_to_hg38(self):
		pass

	def split_by_chrom(self):
		self.vcf_split_paths = {}
		for chrom_name, rows in itt.groupby(self.vcf_ls, key=lambda x: x[0]):
			self.vcf_split_paths[chrom_name] = os.path.join(
				self.vcf_splits_dir,
				'{study_name}_{chrom_name}.vcf'.format(study_name=self.study_name, chrom_name=chrom_name)
				)
			with open(self.vcf_split_paths[chrom_name], 'w') as vcf_file:
				writer = csv.writer(delimiter='\t')
				for row in rows:
					writer.writerow(row)
		return True

	def store(self):
		# store sanitized vcf in db
		pass


class OpenVar:
	def __init__(self, snpeff_path, vcf, annotation='Ensembl+OpenProt'):
		self.snpeff_path  = snpeff_path
		self.snpeff_jar   = os.path.join(snpeff_path, 'snpEff.jar')
		self.snpsift_jar  = os.path.join(snpeff_path, 'SnpSift.jar')

		# make dict annotation -> build
		self.snpeff_build = 'GRCh38.95_refAlt_chr{chrom_name}'

		self.vcf = vcf
		

	def run_snpeff_parallel_pipe(self, nprocs=12):
		pool = multiprocessing.Pool(processes=nprocs)
		pool.map_async(self.run_snpeff, chrom_names)
		return True

	def run_snpeff(self, chrom_name, verbose=False):
		snpEff_chrom_build = self.snpeff_build.format(chrom_name=chrom_name)
		vcf_path           = os.path.join(self.vcf.vcf_split_paths[chrom_name])
		snpEff_logfile     = os.path.join(self.vcf.data_dir, '{}_chr{}_snpEff.log'.format(self.vcf.study_name, chrom_name))

		snpeff_cmd     = self.get_snpeff_cmd(snpEff_chrom_build, vcf_path)
		
		if verbose:
			print('Running SnpEff...')
			print(snpeff_cmd)

		snpeff_subproc = subprocess.Popen(
			snpeff_cmd.split(), 
			shell=False, 
			stdout=subprocess.PIPE,
			stderr=subprocess.STDOUT
			)
		with open(snpEff_logfile, 'w') as f:
			for l in iter(snpeff_subproc.stdout.readline, b''):
				f.write(l.decode())
		snpeff_subproc.wait()

		
		cat_cmd = self.get_cat_cmd()
		perl_cmd = self.get_perl_cmd()
		if verbose:
			print('Formating output...')
			print(cat_cmd)
			print(perl_cmd)
		cat_subproc  = subprocess.Popen(cat_cmd.split(), shell=False, stdout=subprocess.PIPE)
		perl_subproc = subprocess.Popen(perl_cmd.split(), shell=False, stdin=cat_subproc.stdout, stdout=subprocess.PIPE)
		cat_subproc.stdout.close()

		snpsift_cmd = self.get_snpsift_cmd()
		if verbose:
			print(snpsift_cmd)
		annOnePerLine_file = os.path.join(self.vcf.data_dir, '{}_annOnePerLine.tsv'.format(self.vcf.study_name))
		snpsift_subproc    = subprocess.Popen(snpsift_cmd, shell=True, stdin=perl_subproc.stdout, stdout=open(annOnePerLine_file, "w"))

		perl_subproc.stdout.close()
		snpsift_subproc.wait()

		#TODO: catch exceptions

		return True

	def get_snpeff_cmd(self, snpEff_chrom_build, vcf_path):
		cmd = 'java -Xmx12g -jar {snpeff_jar} -v {build} {vcf_path}'.format(
				snpeff_jar = self.snpeff_jar,
				build = snpEff_chrom_build,
				vcf_path = vcf_path
			)
		return cmd

	def get_cat_cmd(self):
		cmd = 'cat {output_vcf}'.format(output_vcf=self.vcf.file_name.replace('.vcf', '.ann.vcf'))
		return cmd

	def get_perl_cmd(self):
		cmd = os.path.join(self.snpeff_path, 'scripts', 'vcfEffOnePerLine.pl')
		return cmd

	def get_snpsift_cmd(self):
		snpsift_fields = 'CHROM POS REF ALT "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].HGVS_P" "ANN[*].FEATUREID" "ANN[*].GENE" "ANN[*].ERRORS" "ANN[*].HGVS_C"'
		cmd ='java -jar {snpsift_jar} extractFields - {snpsift_fields}'.format(
				snpsift_jar = self.snpsift_jar,
				snpsift_fields = snpsift_fields
			)
		return cmd


class OPVReport:
	def __init__(self, opv):
		pass

	def write_tabular(self):
		pass

	def write_html(self):
		pass