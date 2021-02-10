import os
import sys
import csv
import pathlib
import pickle
import pyfaidx
import subprocess
import multiprocessing
import itertools as itt
import numpy as np
from collections import Counter
from matplotlib import pyplot as plt
from pyliftover import LiftOver
import json


maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt / 10)
csv.field_size_limit(sys.maxsize)
gene_len_files = {
    'human': '/open-var-deposit/data/human/gene_lenghts.pkl',
    'mouse': '/open-var-deposit/data/mouse/gene_lenghts.pkl',
    'rat': '/open-var-deposit/data/rat/gene_lenghts.pkl',
    'fruit fly': '/open-var-deposit/data/droso/gene_lenghts.pkl',
}
genome_fastas = {
    'human': '/shared-genomes-folder/human/GRCh38/complete-genome.fa',
    'mouse': '/shared-genomes-folder/human/GRCm38/complete-genome.fa',
    'rat': '',
    'fruit fly': '',
}
chrom_names = [str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']
chrom_set = set(chrom_names)
accepted_bases = {'a', 'c', 'g', 't', 'n', '*'}
vcf_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
impact_levels = {'LOW': 1, 'MODERATE': 2, 'HIGH': 3, 'MODIFIER': 0, 1: 'LOW', 2: 'MODERATE', 3: 'HIGH', 0: 'MODIFIER'}
genome_old_versions = {'hg19': 'hg38', }
annotation_build = {
    ('human', 'OP_Ens'): 'GRCh38.95_refAlt_chr{chrom_name}',
    ('human', 'OP_Ref'): 'GRCh38.p12_chr{chrom_name}',
    ('mouse', 'OP_Ens'): 'GRCm38.95_chr{chrom_name}',
}

class SeqStudy:
    def __init__(self, data_dir, file_name, study_name, results_dir, specie, genome_version, verbose=False):
        self.verbose = verbose
        self.data_dir = mkdir(data_dir)
        self.genome_version = genome_version
        self.specie = specie
        self.results_dir = mkdir(results_dir)
        self.vcf_splits_dir = mkdir(os.path.join(self.results_dir, 'vcf_splits'))
        self.output_dir = mkdir(os.path.join(self.results_dir, 'output'))
        self.file_name = file_name
        self.file_path = os.path.join(data_dir, file_name)
        self.study_name = study_name
        self.warnings = {'unknown chromosomes': [], 'invalid positions': [], 'invalid alleles': []}
        self.file_check = True
        self.parse_vcf()
        print('vcf parsed')
        self.check_vcf_format()
        if self.file_check:
            print('vcf format checked')
            if genome_version in genome_old_versions:
                self.convert_hg19_to_hg38()  # TODO change function to "convert_genome(old_version)"
            self.check_altref_order()
            print('vcf altref allele check')
            self.write_warnings()
            print('vcf warnings written')
            self.split_by_chrom()
            print('vcf chroms splited')

    def parse_vcf(self):
        vcf_ls = []
        with open(self.file_path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for n, row in enumerate(reader):
                if row:
                    vcf_ls.append(row)
        self.vcf_ls = sorted(vcf_ls, key=lambda x: x[0])

    def check_vcf_format(self):
        vcf_ls = []
        for snp in self.vcf_ls:
            snp_line = to_tsv_line(snp)
            snp = dict(zip(vcf_fields, snp))

            # check chromosome names
            if snp['CHROM'].replace('chr', '') not in chrom_set:
                self.warnings['unknown chromosomes'].append(snp['CHROM'])
                continue
            snp['CHROM'] = snp['CHROM'].replace('chr', '')

            # check that position is integer
            try:
                snp['POS'] = int(snp['POS'])
            except:
                self.warnings['invalid positions'].append(snp['POS'])
                continue

            # check alt ref according to VCF 4.2
            if not validate_allele_format(snp['ALT']) and validate_allele_format(snp['REF']):
                self.warnings['invalid alleles'].append(snp_line)
                continue

            vcf_ls.append([snp[field] for field in vcf_fields])

        if not vcf_ls:
            self.file_check = False

        self.vcf_ls = vcf_ls

    def check_altref_order(self):
        vcf_ls = []
        genome = pyfaidx.Fasta(genome_fastas[self.specie], as_raw=True, rebuild=False)
        for snp in self.vcf_ls:
            snp_line = to_tsv_line(snp)
            snp = dict(zip(vcf_fields, snp))
            ref = genome[snp['CHROM']][snp['POS'] - 1]
            if snp['REF'] != ref:
                if snp['ALT'] != ref:
                    self.warnings['invalid alleles'].append(snp_line)
                    continue
                snp['ALT'] = snp['REF']
                snp['REF'] = ref

            vcf_ls.append([snp[field] for field in vcf_fields])
        self.vcf_ls = vcf_ls

    def convert_hg19_to_hg38(self):
        lo_hg38 = LiftOver('hg19', 'hg38')
        vcf_ls = []
        self.warnings['lost at liftOver'] = []
        for snp in self.vcf_ls:
            snp_line = to_tsv_line(snp)
            snp = dict(zip(vcf_fields, snp))
            chrom = 'chr{}'.format(snp['CHROM'])
            lift_hg38 = lo_hg38.convert_coordinate(chrom, snp['POS'])
            if lift_hg38 is not None and lift_hg38:
                hg38_chrom, hg38_pos, strand = lift_hg38[0][0:3]
                snp['CHROM'] = hg38_chrom.replace('chr', '')
                if snp['CHROM'] not in chrom_set:
                    self.warnings['lost at liftOver'].append(snp_line)
                    continue
                snp['POS'] = hg38_pos
            else:
                self.warnings['lost at liftOver'].append(snp_line)
                continue
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

    def write_warnings(self):
        if not any(self.warnings.values()):
            pass
        fpath = os.path.join(self.output_dir, 'warnings.txt')
        with open(fpath, 'w') as f:
            for warning, faults in self.warnings.items():
                if faults:
                    f.write(warning + '\n')
                    for fault in faults:
                        f.write('\t' + fault + '\n')


class OpenVar:
    def __init__(self, snpeff_path, vcf, annotation='OP_Ens'):
        print('creating OPV object...')
        self.snpeff_path = snpeff_path
        self.snpeff_jar = os.path.join(snpeff_path, 'snpEff.jar')
        self.snpsift_jar = os.path.join(snpeff_path, 'SnpSift.jar')
        self.verbose = False
        self.vcf = vcf
        self.specie = vcf.specie
        self.snpeff_build = annotation_build[(self.specie, annotation)]
        self.logs_dir = mkdir(os.path.join(self.vcf.results_dir, 'logs'))
        self.output_dir = self.vcf.output_dir

    def run_snpeff_parallel_pipe(self, nprocs=12):
        pool = multiprocessing.Pool(processes=nprocs)
        r = pool.map(self.run_snpeff, chrom_names)
        pool.close()
        pool.join()

        if all(r):
            return True
        return False

    def run_snpeff(self, chrom_name):
        snpEff_chrom_build = self.snpeff_build.format(chrom_name=chrom_name)
        if chrom_name not in self.vcf.vcf_split_paths:
            if self.verbose:
                print('no variant in chromosome {}'.format(chrom_name))
            return True
        vcf_path = os.path.join(self.vcf.vcf_split_paths[chrom_name])
        vcf_ann_path = vcf_path.replace('.vcf', '.ann.vcf')
        snpEff_logfile = os.path.join(self.logs_dir, '{}_chr{}_snpEff.log'.format(self.vcf.study_name, chrom_name))

        snpeff_cmd = self.get_snpeff_cmd(snpEff_chrom_build, vcf_path)

        if self.verbose:
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
        if self.verbose:
            print('Formating output...')
            print(cat_cmd)
            print(perl_cmd)
        cat_subproc = subprocess.Popen(cat_cmd.split(), shell=False, stdout=subprocess.PIPE)
        perl_subproc = subprocess.Popen(perl_cmd.split(), shell=False, stdin=cat_subproc.stdout, stdout=subprocess.PIPE)
        cat_subproc.stdout.close()

        snpsift_cmd = self.get_snpsift_cmd()
        if self.verbose:
            print(snpsift_cmd)
        annOnePerLine_file = os.path.join(
            self.vcf.vcf_splits_dir, '{study_name}_{chrom_name}_annOnePerLine.tsv'.format(
                study_name=self.vcf.study_name,
                chrom_name=chrom_name
            )
        )
        snpsift_subproc = subprocess.Popen(snpsift_cmd, shell=True, stdin=perl_subproc.stdout,
                                           stdout=open(annOnePerLine_file, "w"))

        perl_subproc.stdout.close()
        snpsift_subproc.wait()

        return True

    def get_snpeff_cmd(self, snpEff_chrom_build, vcf_path):
        cmd = 'java -Xmx12g -jar {snpeff_jar} -v {build} {vcf_path}'.format(
            snpeff_jar=self.snpeff_jar,
            build=snpEff_chrom_build,
            vcf_path=vcf_path
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
        cmd = 'java -jar {snpsift_jar} extractFields - {snpsift_fields}'.format(
            snpsift_jar=self.snpsift_jar,
            snpsift_fields=snpsift_fields
        )
        return cmd


class OPVReport:
    def __init__(self, opv, verbose=False):
        self.opv = opv
        self.verbose = verbose
        self.vcf = opv.vcf
        self.output_dir = self.opv.output_dir
        self.study_name = self.opv.vcf.study_name
        self.specie = opv.specie
        self.annOnePerLine_files = self.list_annOnePerLine_files()
        self.analyzed_variants = self.analyze_all_variants()


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
        with open(annOnePerLine_file_path, 'w') as tab_file:
            writer = csv.writer(tab_file, delimiter='\t')
            for annOnePerLine_file in self.annOnePerLine_files:
                for n, row in enumerate(self.parse_annOnePerLine(annOnePerLine_file, as_dict=True)):
                    if n == 0:
                        writer.writerow(row.keys())
                    writer.writerow(row.values())

        max_impact_file_path = os.path.join(self.output_dir, '{}_max_impact.tsv'.format(self.study_name))
        cols = self.analyzed_variants[0].keys()
        with open(max_impact_file_path, 'w') as tab_file:
            writer = csv.writer(tab_file, delimiter='\t')
            writer.writerow(cols)
            for row in self.analyzed_variants:
                writer.writerow(row.values())

    def dump_to_file(self, var, filename):
        with open(filename, 'w') as f:
            f.write(json.dumps(var, indent=2))

    def compute_summary_stats(self):
        gene_lenghts = pickle.load(open(gene_len_files[self.specie], 'rb'))
        # overall summary
        prot_counts = {
            'alt': len(set(snp['alt_prot_acc'] for snp in self.analyzed_variants if 'IP_' in snp['alt_prot_acc'])),
            'iso': len(set(snp['alt_prot_acc'] for snp in self.analyzed_variants if 'II_' in snp['alt_prot_acc'])),
            'ref': len(set(snp['ref_prot_acc'] for snp in self.analyzed_variants if snp['ref_prot_acc'] != 'null'))
        }

        # chrom level
        snp_set = set([snp['hg38_name'] for snp in self.analyzed_variants])
        snps_per_chroms = Counter([snp.split('_')[0] for snp in snp_set])
        snps_per_chroms = [(chrom, snps_per_chroms[chrom]) for chrom in chrom_names]

        # gene level
        if len(self.analyzed_variants) == len([x for x in self.analyzed_variants if x['gene'] == 'null']):
            gene_snp_rate = []
        else:
            gene_snps_grp = sorted(self.analyzed_variants, key=lambda x: x['gene'])
            #self.dump_to_file(gene_snps_grp, 'var_error.json')
            gene_snp_rate = {gene: len(list(grp)) * 1000 / gene_lenghts[gene] for gene, grp in itt.groupby(gene_snps_grp, key=lambda x: x['gene']) if gene in gene_lenghts}
            gene_snp_rate = sorted(gene_snp_rate.items(), key=lambda x: -x[1])

            fname = os.path.join(self.output_dir, '{}_top_genes_var_rate.svg'.format(self.study_name))
            genes, rates = zip(*gene_snp_rate[:10])  # top ten
            self.generate_bar_chart([genes, rates], 'gene_var_rate', fname)

        # protein level
        if len(self.analyzed_variants) == len([x for x in self.analyzed_variants if x['gene'] == 'null']):
            count_highest= {'alt': 0, 'ref': 0}
            impact_counts = {}
            impact_ann = {}
            fc = {}
        else:
            count_highest = {'alt': sum([snp['alt_max_impact'] > snp['ref_max_impact'] for snp in self.analyzed_variants]), 
                    'ref': sum([snp['ref_max_impact'] > snp['alt_max_impact'] for snp in self.analyzed_variants])}

            impacts = {'ref_all': [x['ref_max_impact'] for x in self.analyzed_variants if x['ref_max_impact'] > -1],
                    'max_all': [max([snp['alt_max_impact'], snp['ref_max_impact']]) for snp in self.analyzed_variants if snp['alt_max_impact'] > -1 or snp['ref_max_impact'] > -1]}

            max_all = dict(zip(range(1, 4), [0]*3))
            max_all.update(dict(Counter(impacts['max_all'])))

            ref_all = dict(zip(range(1, 4), [0]*3))
            ref_all.update(dict(Counter(impacts['ref_all'])))

            impact_ann = {'max_all':max_all, 'ref_all':ref_all}
            fc = {i: max_all[i] / ref_all[i] if ref_all[i] > 0 else 0. for i in range(1, 4)}

            fname = os.path.join(self.output_dir, '{}_impact_foldchange.svg'.format(self.study_name))
            self.generate_bar_chart(fc, 'fold_change', fname)

            impact_counts = dict(zip(range(1, 4), [{'alt': 0, 'ref': 0}, {'alt': 0, 'ref': 0}, {'alt': 0, 'ref': 0}]))
            for snp in self.analyzed_variants:
                if snp['ref_max_impact'] == -1 and snp['alt_max_impact'] == -1:
                    continue
                if snp['ref_max_impact'] and snp['ref_max_impact'] >= snp['alt_max_impact']:
                    impact_counts[snp['ref_max_impact']]['ref'] += 1
                elif snp['alt_max_impact'] and snp['alt_max_impact'] > snp['ref_max_impact']:
                    impact_counts[snp['alt_max_impact']]['alt'] += 1

            ref_impacts = [impact_counts[i]['ref'] for i in range(1, 4)]
            alt_impacts = [impact_counts[i]['alt'] for i in range(1, 4)]
            fname = os.path.join(self.output_dir, '{}_var_per_impact.svg'.format(self.study_name))
            self.generate_bar_chart([ref_impacts, alt_impacts], 'stacked_impact', fname)

        # hotspots on alts
        if len(self.analyzed_variants) == len([x for x in self.analyzed_variants if x['gene'] == 'null']):
            gene_altsnp_rate = {'All nulls': 'No gene consequences for the submitted variants.'}
        else:
            gene_altsnp_rate = {gene: self.count_altsnp_ratio(list(grp)) for gene, grp in itt.groupby(gene_snps_grp, key=lambda x: x['gene']) if gene in gene_lenghts}

        self.summary = {
	        'study_name': self.study_name,
            'Counts summary': {
                'Total number of variants': len(snp_set),
                'Total number of affected genes': len(gene_snp_rate),
                'Total number of affected proteins': sum(prot_counts.values()),
                'Total number of affected reference proteins': prot_counts['ref'],
                'Total number of affected alternative proteins': prot_counts['alt'],
                'Total number of affected novel isoforms': prot_counts['iso'],
            },
            'Chromosome Level': snps_per_chroms,
            'Gene Level': gene_snp_rate,
            'Protein Level': {
                'Number of variants with highest impact on reference proteins': count_highest['ref'],
                'Number of variants with highest impact on alternative proteins': count_highest['alt'],
                'Impact Counts': impact_counts,
                'Impact Annotation':impact_ann,
                'Fold Change': fc,
            },
            'Mutational hotspots on altORFs': gene_altsnp_rate,
        }

        if len(self.analyzed_variants) != len([x for x in self.analyzed_variants if x['gene'] == 'null']):
            fname = os.path.join(self.output_dir, '{}_hotspots_barchart.svg'.format(self.study_name))
            data = zip(*[(gene, counts['ratio_higher_alt'], len(counts['alts'])) for gene, counts in self.summary['Mutational hotspots on altORFs'].items()])
            self.generate_bar_chart(data, 'hotspots_bar', fname)

        fname = os.path.join(self.output_dir, 'summary.pkl')
        pickle.dump(self.summary, open(fname, 'wb'))

    def generate_bar_chart(self, data, chart_type, fname):
        if chart_type == 'gene_var_rate':
            genes, rates = data
            plt.bar(range(1, len(genes) + 1), rates)
            plt.xticks(range(1, len(genes) + 1), genes, rotation='vertical')
            plt.xlabel('Genes')
            plt.ylabel('SNPs per Kb')
            plt.savefig(fname)
            plt.show()

        if chart_type == 'fold_change':
            data = [data[i] for i in range(1, 4)]
            plt.bar(range(1, 4), data)
            plt.xticks(range(1, 4), ['Low', 'Medium', 'High'])
            plt.xlabel('Impact Levels')
            plt.ylabel('Fold Change')
            plt.savefig(fname)
            plt.show()

        if chart_type == 'stacked_impact':
            ref_impacts, alt_impacts = data
            plt.bar(range(1, 4), ref_impacts, label='ref')
            plt.bar(range(1, 4), alt_impacts, bottom=ref_impacts, label='alt')
            plt.xticks(range(1, 4), ['Low', 'Medium', 'High'])
            plt.yscale('log')
            plt.xlabel('Impact Levels')
            plt.ylabel('count SNPs')
            plt.legend()
            plt.savefig(fname)
            plt.show()

        if chart_type == 'hotspots_bar':
            nbin, min_x, max_x = 30, 0., 1.
            genes, freqs, cnt_alts = data
            bins = list(range(1, (nbin + 1)))
            bins = [(min_x + (n - 1) * (max_x / nbin), min_x + n * (max_x / nbin)) for n in bins]
            bin_labels = ['-'.join(['{:.2f}'.format(round(x, 2)) for x in left_right]) for left_right in bins]
            genes_per_bin = {n: [] for n in bins}
            altorf_counts = {n: 0 for n in bins}
            for gene, freq, cnt_alt in zip(genes, freqs, cnt_alts):
                for left, right in bins:
                    if (freq > left) and (freq <= right):
                        genes_per_bin[(left, right)].append(gene)
                        altorf_counts[(left, right)] += cnt_alt

            gene_counts = {n: len(genes_per_bin[n]) for n in bins}
            altorf_per_gene = {n: altorf_counts[n] / gene_counts[n] if gene_counts[n] > 0 else 0. for n in bins}

            Norm = plt.Normalize(min(altorf_per_gene.values()), max(altorf_per_gene.values()))
            val = list(altorf_per_gene.values())
            colors = plt.cm.plasma(Norm(val))

            fig, ax = plt.subplots()
            ax.bar(list(range(1, (nbin + 1), 1)), list(gene_counts.values()), color=colors)
            fig.colorbar(plt.cm.ScalarMappable(norm=Norm, cmap=plt.cm.plasma_r), ax=ax)
            plt.xticks(list(range(1, (nbin + 1))), bin_labels, rotation=90)
            plt.xlabel('binned ratio of SNPs with higher impact in alt')
            plt.ylabel('count genes')
            plt.savefig(fname)
            plt.show()

    def count_altsnp_ratio(self, snp_set):
        cnt_snps = len(set(snp['hg38_name'] for snp in snp_set))
        cnt_alt_snps = len(set(snp['hg38_name'] for snp in snp_set if
                               snp['in_alt'] == 'true' and snp['alt_max_impact'] >= snp['ref_max_impact']))
        alts = list(set(snp['alt_prot_acc'] for snp in snp_set if snp['alt_prot_acc'] != 'null'))
        return {
            'ave_impact': np.mean([x['alt_max_impact'] for x in snp_set]),
            'cnt_snps': cnt_snps,
            'cnt_alt_snps': cnt_alt_snps,
            'ratio_higher_alt': cnt_alt_snps / cnt_snps,
            'alts': alts,
        }

    def analyze_all_variants(self):
        analyzed_variants = []
        for annOnePerLine_file in self.annOnePerLine_files:
            if self.verbose:
                print(annOnePerLine_file)
            snp_effs = self.parse_annOnePerLine(annOnePerLine_file)
            for snp_eff in itt.groupby(sorted(snp_effs, key=lambda x: x[0]), key=lambda x: x[0]):
                analyzed_variants.append(self.analyze_variant(*snp_eff))

        if all([snp['gene'] == 'null' for snp in analyzed_variants]):
            raise Exception('All genes are null!')
        print('All variants analyzed')
        return analyzed_variants

    def list_annOnePerLine_files(self):
        annOnePerLine_files = []
        for f in os.listdir(self.opv.vcf.vcf_splits_dir):
            fpath = os.path.join(self.opv.vcf.vcf_splits_dir, f)
            if os.path.isfile(fpath) and 'annOnePerLine' in fpath:
                annOnePerLine_files.append(fpath)
        return annOnePerLine_files

    def parse_annOnePerLine(self, annOnePerLine_file, as_dict=False):
        fields = ['FEATUREID', 'HGVS_P', 'HGVS_C', 'IMPACT', 'ERRORS', 'GENE']
        with open(annOnePerLine_file, 'r') as f:
            for n, l in enumerate(f):
                ls = l.strip().split('\t')
                if n == 0:
                    keys = ls
                    continue
                line = dict(zip(keys, ls))
                if 'ANN[*].EFFECT' in line:
                    line['ANN[*].EFFECT'] = line['ANN[*].EFFECT'].split('&')
                var_name = '_'.join([line['CHROM'], line['POS'], line['REF'], line['ALT']])
                eff = (var_name, *[line['ANN[*].' + x] if 'ANN[*].' + x in line else 'NA' for x in fields])
                if as_dict:
                    yield line
                else:
                    yield eff

    def analyze_variant(self, variant, effs):
        atts = {
            'hg38_name': variant,
            'in_ref': 'false',
            'in_alt': 'false',
            'ref_max_impact': -1,
            'alt_max_impact': -1,
            'ref_trxpt_acc': 'null',
            'ref_prot_acc': 'null',
            'alt_trxpt_acc': 'null',
            'alt_prot_acc': 'null',
            'ref_hgvs_p': 'null',
            'ref_hgvs_c': 'null',
            'alt_hgvs_p': 'null',
            'alt_hgvs_c': 'null',
            'ref_errs': 'null',
            'alt_errs': 'null',
            'gene': 'null'
        }
        for eff in effs:
            feat_id, hgvs_p, hgvs_c, impact, errs, gene = eff[1:]
            if feat_id_is_ref(feat_id):
                atts['in_ref'] = 'true'
                if impact_levels[impact] > atts['ref_max_impact']:
                    if '@' in feat_id:
                        ref_trxpt_acc, ref_prot_acc = feat_id.split('@')
                        ref_prot_acc = ref_prot_acc.split('.')[0]
                    else:
                        ref_trxpt_acc = feat_id.split('.')[0]
                        ref_prot_acc = ''
                    atts.update({
                        'gene': gene,
                        'ref_trxpt_acc': ref_trxpt_acc,
                        'ref_prot_acc': ref_prot_acc,
                        'ref_hgvs_p': hgvs_p,
                        'ref_hgvs_c': hgvs_c,
                        'ref_errs': errs,
                        'ref_max_impact': impact_levels[impact],
                    })
            elif '^' in feat_id:
                atts['in_alt'] = 'true'
                if impact_levels[impact] > atts['alt_max_impact']:
                    alt_feat_dict = parse_feat_id(feat_id)
                    atts.update(alt_feat_dict)
                    atts.update({
                        'gene': gene,
                        'alt_hgvs_p': hgvs_p,
                        'alt_hgvs_c': hgvs_c,
                        'alt_errs': errs,
                        'alt_max_impact': impact_levels[impact],
                    })
        atts.update({
            'hg38_pos': int(variant.split('_')[1]),
        })
        return atts


def is_synonymous(hgvs_p):
    try:
        hgvs_p = hgvs_p.split('.')[-1]
        return hgvs_p[:3] == hgvs_p[-3:]
    except:
        return False


def parse_feat_id(feat_id):
    if '^' in feat_id:
        trxpt, acc = feat_id.split('^')
    if acc.count('_')>1:
        acc = '_'.join(acc.split('_')[:2])
    return {'alt_trxpt_acc': trxpt, 'alt_prot_acc': acc}


def feat_id_is_ref(feat_id):
    if 'IP_' not in feat_id and '^' not in feat_id:
        return True
    return False


def mkdir(path):
    try:
        pathlib.Path(path).mkdir(parents=True, exist_ok=True)
    except:
        print("can't create dir: {}".format(path))
    return path


def to_tsv_line(ls):
    return '\t'.join([str(x) for x in ls])


def validate_allele_format(allele):
    for a in allele.split(','):
        for nt in a:
            if nt.lower() not in accepted_bases:
                return False
    return True
