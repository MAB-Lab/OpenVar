import os
import sys
import csv
import json
import pathlib
import pickle
import pyfaidx
import subprocess
import multiprocessing
import itertools as itt
from collections import Counter
from matplotlib import pyplot as plt
from pyliftover import LiftOver

maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)
csv.field_size_limit(sys.maxsize)

hg38_genome = pyfaidx.Fasta('/home/xroucou_group/genomes/human/GRCh38/complete-genome.fa', as_raw=True)
chrom_names = [str(x) for x in range(1,23)] + ['X', 'Y', 'MT']
vcf_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
impact_levels = {'LOW':1, 'MODERATE':2, 'HIGH':3, 'MODIFIER':0, 1:'LOW', 2:'MODERATE', 3:'HIGH', 0:'MODIFIER'}
prot_gene_dict = pickle.load(open('OP1.6_prot_gene_dict.pkl', 'rb'))
gene_lenghts = pickle.load(open('gene_lenghts.pkl', 'rb'))
genome_old_versions = {'hg19':'hg38', 'mm39':'mm10', 'rn6':'rn5', 'dm6':'dm5'}

class SeqStudy:
	def __init__(self, data_dir, file_name, study_name, results_dir, genome_version):
		self.data_dir    = mkdir(data_dir)
		self.results_dir = mkdir(results_dir)
		self.vcf_splits_dir = mkdir(os.path.join(self.results_dir, 'vcf_splits'))
		self.output_dir = mkdir(os.path.join(self.results_dir, 'output'))
		self.file_name  = file_name
		self.file_path  = os.path.join(data_dir, file_name)
		self.study_name = study_name
		self.warnings   = []
		self.file_check = True # TODO switch flag at proper check failures...
		self.parse_vcf()
		self.check_vcf_format()
		if genome_version in genome_old_versions:
			self.convert_hg19_to_hg38() #TODO change funtion to convert_genome(old_version)
		self.check_altref_order()
		self.write_warnings()
		self.split_by_chrom()

	def parse_vcf(self):
		# accept comma delimiters
		vcf_ls = []
		with open(self.file_path, 'r') as f:
			reader = csv.reader(f, delimiter='\t')
			for n, row in enumerate(reader):
				vcf_ls.append(row)
		self.vcf_ls = sorted(vcf_ls, key=lambda x: x[0])

	def check_vcf_format(self):
		chrom_set = set(chrom_names)
		vcf_ls = []
		# is any chrom in chrome set?
		# for these, is pos an integer?
		# check alt ref accoring to VCF 4.2
		for snp in self.vcf_ls:
			snp_line = to_tsv_line(snp)
			snp = dict(zip(vcf_fields, snp))

			# check chrom names
			snp['CHROM'] = snp['CHROM'].replace('chr', '')
			if snp['CHROM'] not in chrom_set:
				self.warnings.append('Unknown chromosome: {}'.format(snp['CHROM']))
				continue

			# check that position is integer
			try:
				snp['POS'] = int(snp['POS'])
			except:
				self.warnings.append('position is not an integer: {}'.format(snp_line))
				continue
			vcf_ls.append([snp[field] for field in vcf_fields])
		self.vcf_ls = vcf_ls

	def check_altref_order(self):
		vcf_ls = []
		for snp in self.vcf_ls:
			snp_line = to_tsv_line(snp)
			snp = dict(zip(vcf_fields, snp))
			ref = hg38_genome[snp['CHROM']][snp['POS']-1]
			if snp['REF'] != ref:
				if snp['ALT'] == ref
					snp['REF'] = snp['ALT']
					snp['ALT'] = ref
					self.warnings.append('switched REF and ALT alleles to match ref genome: {}'.format(snp_line))
			vcf_ls.append([snp[field] for field in vcf_fields])
		self.vcf_ls = vcf_ls
		# if not any match, raise error

	def convert_hg19_to_hg38(self):
		lo_hg38 = LiftOver('hg19', 'hg38')
		vcf_ls = []
		for snp in self.vcf_ls:
			snp_line = to_tsv_line(snp)
			snp = dict(zip(vcf_fields, snp))
			chrom = 'chr{}'.format(snp['CHROM'])
			lift_hg38 = lo_hg38.convert_coordinate(chrom, snp['POS'])
			if lift_hg38 is not None and lift_hg38:
				hg38_chrom, hg38_pos, strand = lift_hg38[0][0:3]
				snp['CHROM']  = hg38_chrom
				snp['POS']  = hg38_pos
			else:
				self.warnings.append('lost at liftOver conversion: {}'.format(snp_line))
			vcf_ls.append([snp[field] for field in vcf_fields])
		self.vcf_ls = vcf_ls

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

	def write_warnings(self):
		fpath = os.path.join(self.output_dir, 'warnings.txt') 
		with open(fpath, 'w') as f:
			for warning in self.warnings:
				f.write(warning+'\n')

class OpenVar:
	def __init__(self, snpeff_path, vcf, annotation='OP_Ens'):
		self.snpeff_path  = snpeff_path
		self.snpeff_jar   = os.path.join(snpeff_path, 'snpEff.jar')
		self.snpsift_jar  = os.path.join(snpeff_path, 'SnpSift.jar')

		# make dict annotation -> build
		self.snpeff_build = 'GRCh38.95_refAlt_chr{chrom_name}'
		self.vcf = vcf
		
		self.logs_dir = mkdir(os.path.join(self.vcf.results_dir, 'logs'))
		self.output_dir = self.vcf.output_dir

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
		snpEff_logfile     = os.path.join(self.logs_dir, '{}_chr{}_snpEff.log'.format(self.vcf.study_name, chrom_name))

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
		self.output_dir = self.opv.output_dir
		self.study_name = self.opv.vcf.study_name
		self.parse_annOnePerLine()

	def aggregate_annotated_vcf(self):
		split_ann_vcfs = []
		for f in os.listdir(self.opv.vcf.vcf_splits_dir):
		    fpath = os.path.join(self.opv.vcf.vcf_splits_dir, f)
		    if os.path.isfile(fpath) and '.ann.vcf' in fpath:
		        split_ann_vcfs.append(fpath)

		ann_vcf_file_path = os.path.join(self.output_dir, self.vcf.file_name.replace('.vcf', '.ann.vcf'))        
		with open(ann_vcf_file_path, 'w') as ann_vcf_file:
			writer = csv.writer(ann_vcf_file, delimiter='\t')
			for split_ann_vcf_file in split_ann_vcfs:
				with open(split_ann_vcf_file, 'r') as fread:
					reader = csv.reader(fread, delimiter='\t')
					for row in reader:
						writer.writerow(row)

	def write_tabular(self):
		annOnePerLine_file_path = os.path.join(self.output_dir, '{}_annOnePerLine.tsv'.format(self.study_name))
		cols = self.annOnePerLine[0].keys()
		with open(annOnePerLine_file_path, 'w') as tab_file:
			writer = csv.writer(tab_file, delimiter='\t')
			writer.writerow(cols)
			for row in self.annOnePerLine:
				writer.writerow(row.values())

	def compute_summary_stats(self):
		# overall summary
		prot_counts = {
			'alt' :len(set(snp['alt_prot_acc'] for snp in self.analyzed_variants if 'IP_' in snp['alt_prot_acc'])),
			'iso' :len(set(snp['alt_prot_acc'] for snp in self.analyzed_variants if 'II_' in snp['alt_prot_acc'])),
			'ref' :len(set(snp['ref_prot_acc'] for snp in self.analyzed_variants))
		}

		# chrom level
		snp_set = set([snp['hg38_name'] for snp in self.analyzed_variants])
		snps_per_chroms = Counter([snp.split('_')[0] for snp in snp_set])
		snps_per_chroms = [(chrom, snps_per_chroms[chrom]) for chrom in chrom_names]

		# gene level
		gene_snps_grp = sorted(self.analyzed_variants, key=lambda x: x['gene'])
		gene_snp_rate = {gene:len(list(grp))*1000/gene_lenghts[gene] for gene, grp in itt.groupby(gene_snps_grp, key=lambda x: x['gene']) if gene in gene_lenghts}
		gene_snp_rate = sorted(gene_snp_rate.items(), key=lambda x: -x[1])

		fname = os.path.join(self.output_dir, '{}_top_genes_var_rate.svg'.format(self.study_name))
		genes, rates = zip(*gene_snp_rate[:10]) # top ten
		self.generate_bar_chart([genes, rates], 'gene_var_rate', fname)

		# protein level
		count_higest = {
			'alt':sum([snp['alt_max_impact']>snp['ref_max_impact'] for snp in self.analyzed_variants]),
			'ref':sum([snp['ref_max_impact']>snp['alt_max_impact'] for snp in self.analyzed_variants])
		}

		impacts = {
			'ref_all':[x['ref_max_impact'] for x in self.analyzed_variants if x['ref_max_impact']>-1],
			'max_all':[max([snp['alt_max_impact'], snp['ref_max_impact']]) for snp in self.analyzed_variants if snp['alt_max_impact']>-1 or snp['ref_max_impact']>-1],
		}

		max_all = dict(Counter(impacts['max_all']))
		ref_all = dict(Counter(impacts['ref_all']))
		fc = {i:max_all[i]/ref_all[i] for i in range(1,4)}
		fname = os.path.join(self.output_dir, '{}_impact_foldchange.svg'.format(self.study_name))
		self.generate_bar_chart(fc, 'fold_change', fname)


		impact_counts = dict(zip(range(1,4), [{'alt':0, 'ref':0}, {'alt':0, 'ref':0}, {'alt':0, 'ref':0}]))
		for snp in self.analyzed_variants:
			if snp['ref_max_impact']==-1 and snp['alt_max_impact']==-1: continue
			if snp['ref_max_impact']>=snp['alt_max_impact']:
				impact_counts[snp['ref_max_impact']]['ref'] += 1
			elif snp['alt_max_impact']>snp['ref_max_impact']:
				impact_counts[snp['alt_max_impact']]['alt'] += 1

		ref_impacts = [impact_counts[i]['ref'] for i in range(1,4)]
		alt_impacts = [impact_counts[i]['alt'] for i in range(1,4)]
		fname = os.path.join(self.output_dir, '{}_var_per_impact.svg'.format(self.study_name))
		self.generate_bar_chart([ref_impacts, alt_impacts], 'stacked_impact', fname)


		# hotspots on alts
		gene_altsnp_rate = {gene:self.count_altsnp_ratio(list(grp)) for gene, grp in itt.groupby(gene_snps_grp, key=lambda x: x['gene']) if gene in gene_lenghts}

		self.summary = {
			'Counts summary':{
				'Total number of variants':len(snp_set),
				'Total number of affected genes':len(gene_snp_rate),
				'Total number of affected proteins':sum(prot_counts.values()),
				'Total number of affected reference proteins':prot_counts['ref'],
				'Total number of affected alternative proteins':prot_counts['ref'],
				'Total number of affected novel isoforms':prot_counts['iso'],
			},
			'Chromosome Level':snps_per_chroms,
			'Gene Level':gene_snp_rate,
			'Protein Level':{
				'Number of variants with highest impact on reference proteins':count_higest['ref'],
				'Number of variants with highest impact on alternative proteins':count_higest['alt'],
				'Impact Counts': impact_counts,
				'Fold Change':fc,
			},
			'Mutational hotspots on altORFs':gene_altsnp_rate,
		}

		fname = os.path.join(self.output_dir, '{}_hotspots_barchart.svg'.format(self.study_name))
		data = zip(*[(gene, counts['ratio_higher_alt'], len(counts['alts'])) for gene, counts in self.summary['Mutational hotspots on altORFs'].items()])
		self.generate_bar_chart(data, 'hotspots_bar', fname)

		fname = os.path.join(self.output_dir, '{}_summary.pkl'.format(self.study_name))
		pickle.dump(self.summary, open(fname, 'rb'))

	def generate_bar_chart(self, data, chart_type, fname):
		if chart_type=='gene_var_rate':
			genes, rates = data
			plt.bar(range(1,11), rates)
			plt.xticks(range(1,11), genes, rotation='vertical')
			plt.xlabel('Genes')
			plt.ylabel('SNPs per Kb')
			plt.savefig(fname)
			plt.show()

		if chart_type=='fold_change':
			data = [data[i] for i in range(1,4)]
			plt.bar(range(1,4), data)
			plt.xticks(range(1,4), ['Low', 'Medium', 'High'])
			plt.xlabel('Impact Levels')
			plt.ylabel('Fold Change')
			plt.savefig(fname)
			plt.show()

		if chart_type=='stacked_impact':
			ref_impacts, alt_impacts = data
			plt.bar(range(1,4), ref_impacts)
			plt.bar(range(1,4), alt_impacts, bottom=ref_impacts)
			plt.xticks(range(1,4), ['Low', 'Medium', 'High'])
			plt.yscale('log')
			plt.xlabel('Impact Levels')
			plt.ylabel('count SNPs')
			plt.legend()
			plt.savefig(fname)
			plt.show()

		if chart_type=='hotspots_bar':
			nbin = 30
			genes, freqs, cnt_alts = data
			bins = list(range(1, (nbin + 1)))
			bins = [(min_x + (n-1)*(max_x/nbin), min_x + n*(max_x/nbin)) for n in bins]
			bin_labels = ['-'.join(['{:.2f}'.format(round(x,2)) for x in left_right]) for left_right in bins]
			genes_per_bin = {n:[] for n in bins}
			altorf_counts = {n:0 for n in bins}
			for gene, freq, cnt_alt in zip(genes, freqs, cnt_alts):
				for left, right in bins:
					if (freq > left) and (freq <= right):
						genes_per_bin[(left, right)].append(gene)
						altorf_counts[(left, right)] += cnt_alt
			            
			gene_counts = {n:len(genes_per_bin[n]) for n in bins}
			altorf_per_gene = {n:altorf_counts[n]/gene_counts[n] for n in bins}

			Norm = plt.Normalize(min(altorf_per_gene.values()), max(altorf_per_gene.values()))
			val = list(altorf_per_gene.values())
			colors = plt.cm.plasma(Norm(val))

			fig, ax = plt.subplots()
			bar = ax.bar(list(range(1, (nbin + 1), 1)), list(gene_counts.values()), color = colors)
			fig.colorbar(plt.cm.ScalarMappable(norm = Norm, cmap = plt.cm.plasma), ax = ax)
			plt.xticks(list(range(1, (nbin + 1))), bin_labels, rotation=90)
			plt.xlabel('binned ratio of SNPs with higher impact in alt')
			plt.ylabel('count genes')
			plt.savefig(fname)
			plt.show()

	def count_altsnp_ratio(self, snp_set):
		cnt_snps = len(set(snp['hg38_name'] for snp in snp_set))
		cnt_alt_snps = len(set(snp['hg38_name'] for snp in snp_set if snp['in_alt']=='true' and snp['alt_max_impact']>snp['ref_max_impact']))
		alts = list(set(snp['alt_prot_acc'] for snp in snp_set if snp['alt_prot_acc']!='null'))
		return {
			#'snp_set': snp_set,
			'cnt_snps':cnt_snps,
			'cnt_alt_snps':cnt_alt_snps,
			'score':(cnt_alt_snps/cnt_snps)*cnt_alt_snps,
			'ratio_higher_alt':cnt_alt_snps/cnt_snps,
			'alts':alts
		}

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
			if hgvs_p:
				if 'ENST' in feat_id and '@' in feat_id:
					atts['in_ref'] = 'true'
					if impact_levels[impact] and impact_levels[impact] > atts['ref_max_impact']:
						ref_trxpt_acc, ref_prot_acc = feat_id.split('@')
						ref_prot_acc = ref_prot_acc.split('.')[0]
						gene = prot_gene_dict[ref_prot_acc]
						atts.update({
							'gene'           : gene,
							'ref_trxpt_acc'  : ref_trxpt_acc,
							'ref_prot_acc'   : ref_prot_acc,
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
	if '^' in feat_id:
		trxpt, acc = feat_id.split('^')
	return {'alt_trxpt_acc':trxpt, 'alt_prot_acc':acc, 'gene': prot_gene_dict[acc]}

def mkdir(path):
	try:
		pathlib.Path(path).mkdir(parents=True, exist_ok=True)
	except:
		print("can't create dir: {}".format(path))
	return path

def to_tsv_line(ls):
	return '\t'.join([str(x) for x in ls])