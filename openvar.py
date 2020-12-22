import os
import csv
import subprocess
import multiprocessing
import itertools as itt

chrom_names = set([str(x) for x in range(1,23)] + ['X', 'Y', 'MT'])
impact_levels = {'LOW':1, 'MODERATE':2, 'HIGH':3, 'MODIFIER':0, 1:'LOW', 2:'MODERATE', 3:'HIGH', 0:'MODIFIER'}


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
			for n, row in enumerate(reader):
				vcf_ls.append(row)
		self.vcf_ls = sorted(vcf_ls, key=lambda x: x[0])

	def run_all_checks(self):
		# check ref-alt allele order (GNOMAD)
		#for snp in self.vcf_ls:

		# check genome version - convert?
		
		# check chrome names

		# return proper feedback in case of non-compliance
		pass

	def check_alt_ref_allele(self):

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
				writer = csv.writer(vcf_file, delimiter='\t')
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
		r = pool.map(self.run_snpeff, chrom_names)
		pool.close()
		pool.terminate()
		pool.join()
		if all(r):
			return True
		return False

	def run_snpeff(self, chrom_name, verbose=False):
		snpEff_chrom_build = self.snpeff_build.format(chrom_name=chrom_name)
		if chrom_name not in self.vcf.vcf_split_paths:
			print('no variant in chromosome {}'.format(chrom_name))
			return True
		vcf_path           = os.path.join(self.vcf.vcf_split_paths[chrom_name])
		vcf_ann_path       = vcf_path.replace('.vcf', '.ann.vcf')
		snpEff_logfile     = os.path.join(self.vcf.data_dir, '{}_chr{}_snpEff.log'.format(self.vcf.study_name, chrom_name))

		snpeff_cmd     = self.get_snpeff_cmd(snpEff_chrom_build, vcf_path)
		
		if verbose:
			print('Running SnpEff...')
			print(snpeff_cmd)

		snpeff_subproc = subprocess.Popen(
			snpeff_cmd.split(), 
			shell=False, 
			stdout=subprocess.PIPE,
			stderr=open(snpEff_logfile, 'w')
			)
		with open(vcf_ann_path, 'w') as f:
			for l in iter(snpeff_subproc.stdout.readline, b''):
				f.write(l.decode())
		snpeff_subproc.wait()

		
		cat_cmd = self.get_cat_cmd(vcf_ann_path)
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
		annOnePerLine_file = os.path.join(
			self.vcf.vcf_splits_dir, '{study_name}_{chrom_name}_annOnePerLine.tsv'.format(
				study_name=self.vcf.study_name,
				chrom_name=chrom_name
				)
			)
		snpsift_subproc    = subprocess.Popen(snpsift_cmd, shell=True, stdin=perl_subproc.stdout, stdout=open(annOnePerLine_file, "w"))

		perl_subproc.stdout.close()
		snpsift_subproc.wait()

		return True

	def get_snpeff_cmd(self, snpEff_chrom_build, vcf_path):
		cmd = 'java -Xmx12g -jar {snpeff_jar} -v {build} {vcf_path}'.format(
				snpeff_jar = self.snpeff_jar,
				build = snpEff_chrom_build,
				vcf_path = vcf_path
			)
		return cmd

	def get_cat_cmd(self, vcf_ann_path):
		cmd = 'cat {output_vcf}'.format(output_vcf=vcf_ann_path)
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
		self.opv = opv
		self.vcf = opv.vcf
		self.data_dir = self.opv.vcf.data_dir
		self.study_name = self.opv.vcf.study_name
		self.parse_annOnePerLine()

	def aggregate_annotated_vcf(self):
		split_ann_vcfs = []
		for f in os.listdir(self.opv.vcf.vcf_splits_dir):
		    fpath = os.path.join(self.opv.vcf.vcf_splits_dir, f)
		    if os.path.isfile(fpath) and '.ann.vcf' in fpath:
		        split_ann_vcfs.append(fpath)

		ann_vcf_file_path = os.path.join(self.data_dir, self.file_name.replace('.vcf', '.ann.vcf'))        
		with open(ann_vcf_file_path, 'w') as ann_vcf_file:
			writer = csv.writer(ann_vcf_file, delimiter='\t')
			for split_ann_vcf_file in split_ann_vcfs:
				with open(split_ann_vcf_file, 'r') as fread:
					reader = csv.reader(fread, delimiter='\t')
					for row in reader:
						writer.writerow(row)

	def write_tabular(self):
		annOnePerLine_file_path = os.path.join(self.data_dir, '{}_annOnePerLine.tsv'.format(self.study_name))
		cols = self.annOnePerLine[0].keys()
		with open(annOnePerLine_file_path, 'w') as tab_file:
			writer = csv.writer(tab_file, delimiter='\t')
			witer.writerow(cols)
			for row in self.annOnePerLine:
				witer.writerow(row.values())

	def write_json_report(self):
		pass

	def analyze_all_variants(self):
		snps = []
		for snp in self.annOnePerLine:
			var_name = '_'.join([snp['CHROM'], snp['POS'], snp['REF'], snp['ALT']])
			snps.append(
				(var_name, snp['ANN[*].FEATUREID'], snp['ANN[*].HGVS_P'], snp['ANN[*].HGVS_C'], snp['ANN[*].IMPACT'], snp['ANN[*].ERRORS'])
				)
		analyzed_variants = []
		for snp in itt.groupby(snps, key=lambda x: x[0]):
			analyzed_variants.append(self.analyze_variant(*snp))
		self.analyzed_variants = analyzed_variants

	def parse_annOnePerLine(self):
		annOnePerLines = []
		for f in os.listdir(self.opv.vcf.vcf_splits_dir):
		    fpath = os.path.join(self.opv.vcf.vcf_splits_dir, f)
		    if os.path.isfile(fpath) and 'annOnePerLine' in fpath:
		        annOnePerLines.append(fpath)
		lines = []
		for annOnePerLine in annOnePerLines:
		    with open(annOnePerLine, 'r') as f:
		        for n,l in enumerate(f):
		            ls = l.strip().split('\t')
		            if n==0:
		                keys=ls
		                continue
		            line=dict(zip(keys, ls))
		            line['ANN[*].EFFECT'] = line['ANN[*].EFFECT'].split('&')
		            lines.append(line)
		self.annOnePerLine = lines

	def analyze_variant(self, variant, effs, debug=False):
	    atts = {
	        'hg38_name'      : variant,
	        'in_ref'         : 'false',
	        'in_alt'         : 'false',
	        'ref_max_impact' : -1,
	        'alt_max_impact' : -1,
	        'ref_trxpt_acc'  : 'null',
	        'ref_prot_acc'   : 'null',
	        'alt_trxpt_acc'  : 'null',
	        'alt_prot_acc'   : 'null',
	        'ref_hgvs_p'     : 'null',
	        'ref_hgvs_c'     : 'null',
	        'alt_hgvs_p'     : 'null',
	        'alt_hgvs_c'     : 'null',
	        'ref_errs'       : 'null',
	        'alt_errs'       : 'null',
	        'gene'           : 'null'
	           }
	    for eff in effs:
	        feat_id, hgvs_p, hgvs_c, impact, errs = eff[1:]
	        if debug:
	            print('-----------')
	            print(feat_id, hgvs_p, hgvs_c, impact, errs)
	            print("hgvs_p: "+hgvs_p, "'^' in feat_id: "+str('^' in feat_id),"max_impact: "+str(impact_levels[impact] > atts['alt_max_impact']))
	        if hgvs_p:
	            if 'ENST' in feat_id and '@' in feat_id:
	                atts['in_ref'] = 'true'
	                if impact_levels[impact] and impact_levels[impact] > atts['ref_max_impact']:
	                    atts.update({
	                        'ref_trxpt_acc'  : feat_id,
	                        'ref_hgvs_p'     : hgvs_p,
	                        'ref_hgvs_c'     : hgvs_c,
	                        'ref_errs'       : errs,
	                        'ref_max_impact' : impact_levels[impact]
	                    })
	            elif '^' in feat_id:
	                atts['in_alt'] = 'true'
	                if impact_levels[impact] and impact_levels[impact] > atts['alt_max_impact']:
	                    alt_feat_dict = parse_feat_id(feat_id)
	                    atts.update(alt_feat_dict)
	                    atts.update({
	                        'alt_hgvs_p'     : hgvs_p,
	                        'alt_hgvs_c'     : hgvs_c,
	                        'alt_errs'       : errs,
	                        'alt_max_impact' : impact_levels[impact],
	                    })
	    atts.update({
	        'hg38_pos'        : int(variant.split('_')[1]),
	    })
	    return atts

def is_synonymous(hgvs_p):
	try:
	    hgvs_p = hgvs_p.split('.')[-1]
	    return hgvs_p[:3] == hgvs_p[-3:]
	except:
	        return False

def parse_feat_id(feat_id, gene_dict=None):
	trxpt, acc = feat_id.split('^')
	return {'alt_trxpt_acc':trxpt, 'alt_prot_acc':acc, 'gene': op_prot_gene[acc]}
