import os
import sys
import csv
import json
import pathlib
import pickle
import pyfaidx
import subprocess
import threading
import numpy as np
import itertools as itt
from shutil import copyfile
from collections import Counter
from matplotlib import pyplot as plt
from pyliftover import LiftOver

maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt / 10)
csv.field_size_limit(sys.maxsize)
gene_len_files = {
    ('human', 'Ensembl'): '/open-var-deposit/data/human/gene_lenghts_ensembl.pkl',
    ('mouse', 'Ensembl'): '/open-var-deposit/data/mouse/gene_lenghts_ensembl.pkl',
    ('human', 'RefSeq'): '/open-var-deposit/data/human/gene_lenghts_refseq.pkl',
    ('mouse', 'RefSeq'): '/open-var-deposit/data/mouse/gene_lenghts_refseq.pkl',
}
genome_fastas = {
    'human': '/shared-genomes-folder/human/GRCh38/complete-genome.fa',
    'mouse': '/shared-genomes-folder/mouse/GRCm38/complete-genome.fa',
    'rat': '',
    'fruit fly': '',
}
chrom_names = {
    'human': [str(x) for x in range(1, 23)] + ['X', 'Y', 'MT'],
    'mouse': [str(x) for x in range(1, 20)] + ['X', 'Y', 'MT']
}
chrom_set = {sp:set(chrom_names[sp]) for sp in chrom_names}

accepted_bases = {'a', 'c', 'g', 't', 'n', '*'}
vcf_fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT']
impact_levels = {'LOW': 1, 'MODERATE': 2, 'HIGH': 3, 'MODIFIER': 0, 1: 'LOW', 2: 'MODERATE', 3: 'HIGH', 0: 'MODIFIER'}
genome_old_versions = {'hg19': 'hg38', }
annotation_build = {
    ('human', 'OP_Ensembl'): 'GRCh38.95_refAlt_chr{chrom_name}',
    ('human', 'OP_RefSeq'): 'GRCh38.p12_chr{chrom_name}',
    ('human', 'Ensembl'): 'GRCh38.95',
    ('human', 'RefSeq'): 'GRCh38.p12',
    ('mouse', 'OP_Ensembl'): 'GRCm38.95_chr{chrom_name}',
    ('mouse', 'OP_RefSeq'): 'GRCm38.p6_chr{chrom_name}',
    ('mouse', 'Ensembl'): 'GRCm38.95',
    ('mouse', 'RefSeq'): 'GRCm38.p6',
}

class SeqStudy:
    def __init__(self, data_dir, file_name, study_name, results_dir, specie, genome_version, annotation, verbose=False):
        self.verbose = verbose
        self.specie = specie
        self.genome_version = genome_version
        self.annotation = annotation
        self.study_name = study_name
        self.file_name = file_name
        self.data_dir = mkdir(data_dir)
        self.results_dir = mkdir(results_dir)
        self.output_dir = mkdir(os.path.join(self.results_dir, 'output'))
        self.file_path = os.path.join(data_dir, file_name)
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
            if 'OP' in self.annotation: # split vcf by chrom if using OP annoations
                self.vcf_splits_dir = mkdir(os.path.join(self.results_dir, 'vcf_splits'))
                self.split_by_chrom()
                print('vcf chroms splited')
            else:
                self.single_vcf_path = os.path.join(self.results_dir, self.study_name + '.vcf')
                self.write_vcf(self.single_vcf_path, self.vcf_ls)

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
            if snp['CHROM'].replace('chr', '') not in chrom_set[self.specie]:
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
            ref_alt = ref
            if (',' in snp['ALT']) or (',' in snp['REF']):
                self.warnings['invalid alleles'].append(snp_line)
                continue
            if len(snp['REF']) > 1:
                pos_end_r = snp['POS'] - 1 + len(snp['REF'])
                ref = genome[snp['CHROM']][snp['POS'] - 1:pos_end_r]
            if len(snp['ALT']) > 1:
                pos_end_a = snp['POS'] - 1 + len(snp['ALT'])
                ref_alt = genome[snp['CHROM']][snp['POS'] - 1:pos_end_a]
            if snp['REF'] != ref:
                if snp['ALT'] != ref_alt:
                    self.warnings['invalid alleles'].append(snp_line)
                    continue
                snp['ALT'] = snp['REF']
                snp['REF'] = ref_alt

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
                if snp['CHROM'] not in chrom_set[self.specie]:
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
            self.write_vcf(self.vcf_split_paths[chrom_name], rows)

    def write_vcf(self, fpath, vcf_rows):
        with open(fpath, 'w') as vcf_file:
            writer = csv.writer(vcf_file, delimiter='\t')
            for row in vcf_rows:
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


class Worker:
    def __init__(self, function, args):
        self.result = None
        self.args = args

        def work():
            self.result = function(self.args)

        self.thread = threading.Thread(target = work)

class OpenVar:
    def __init__(self, snpeff_path, vcf):
        self.vcf = vcf
        self.specie = vcf.specie
        self.annotation = vcf.annotation
        self.snpeff_build = annotation_build[(self.specie, self.annotation)]
        self.snpeff_path = snpeff_path
        self.snpeff_jar = os.path.join(snpeff_path, 'snpEff.jar')
        self.snpsift_jar = os.path.join(snpeff_path, 'SnpSift.jar')
        self.verbose = False
        self.logs_dir = mkdir(os.path.join(self.vcf.results_dir, 'logs'))
        self.output_dir = self.vcf.output_dir

    def run_snpeff_parallel_pipe(self, nprocs=12):
        workers = [ Worker(self.run_snpeff_chromosome, chrom_name) for chrom_name in chrom_names[self.specie] ]
        
        for w in workers:
            w.thread.start()

        for w in workers:
            w.thread.join()

        for w in workers:
            if w.result is None or w.result is False:
                return False

        return True

    def run_snpeff_chromosome(self, chrom_name):
        snpEff_chrom_build = self.snpeff_build.format(chrom_name=chrom_name)
        if chrom_name not in self.vcf.vcf_split_paths:
            if self.verbose:
                print('no variant in chromosome {}'.format(chrom_name))
            return True
        vcf_path = os.path.join(self.vcf.vcf_split_paths[chrom_name])
        return self.run_snpeff(vcf_path, snpEff_chrom_build)

    def run_snpeff(self, vcf_path, build):
        if build in {'Ensembl', 'RefSeq'}:
            build = annotation_build[(self.specie, build)]
        vcf_ann_path = vcf_path.replace('.vcf', '.ann.vcf')
        snpEff_logfile = os.path.join(self.logs_dir, '{}_{}_snpEff.log'.format(self.vcf.study_name, build))

        snpeff_cmd = self.get_snpeff_cmd(build, vcf_path)

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
        if hasattr(self.vcf, 'vcf_splits_dir'):
            vcf_dir = self.vcf.vcf_splits_dir
        else:
            vcf_dir = self.vcf.results_dir
        annOnePerLine_file = os.path.join(
            vcf_dir, '{study_name}_{build}_annOnePerLine.tsv'.format(
                study_name=self.vcf.study_name,
                build=build
            )
        )
        snpsift_subproc = subprocess.Popen(snpsift_cmd, shell=True, stdin=perl_subproc.stdout,
                                           stdout=open(annOnePerLine_file, "w"))

        perl_subproc.stdout.close()
        snpsift_subproc.wait()

        return True

    def get_snpeff_cmd(self, build, vcf_path):
        cmd = 'java -Xmx12g -jar {snpeff_jar} -v {build} {vcf_path}'.format(
            snpeff_jar=self.snpeff_jar,
            build=build,
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
        if 'OP_' not in self.opv.annotation:
            single_vcf_ann = self.opv.vcf.single_vcf_path.replace('.vcf', '.ann.vcf')
            dst = os.path.join(self.output_dir, os.path.split(single_vcf_ann)[-1])
            copyfile(single_vcf_ann, dst)
            return True

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
        if 'OP_' in self.opv.annotation:
            with open(annOnePerLine_file_path, 'w') as tab_file:
                writer = csv.writer(tab_file, delimiter='\t')
                header_written = False
                for annOnePerLine_file in self.annOnePerLine_files:
                    for n, row in enumerate(self.parse_annOnePerLine(annOnePerLine_file, as_dict=True)):
                        if n == 0 and not header_written:
                            writer.writerow(row.keys())
                            header_written = True
                        writer.writerow(row.values())
        else:
            copyfile(self.annOnePerLine_files[0], annOnePerLine_file_path)

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

    def compute_chrom_gene_level_stats(self, write_summary_pkl=False):
        gene_lenghts = pickle.load(open(gene_len_files[(self.specie, self.opv.annotation.replace('OP_', ''))], 'rb'))

        snp_set = set([snp['hg38_name'] for snp in self.analyzed_variants])
        snps_per_chroms = Counter([snp.split('_')[0] for snp in snp_set])
        snps_per_chroms = [(chrom, snps_per_chroms[chrom]) for chrom in chrom_names[self.specie]]

        gene_snps_grp = sorted(self.analyzed_variants, key=lambda x: x['gene'])
        gene_snp_rate = {gene: len(list(grp)) * 1000 / gene_lenghts[gene] for gene, grp in itt.groupby(gene_snps_grp, key=lambda x: x['gene']) if gene in gene_lenghts}
        gene_snp_rate = sorted(gene_snp_rate.items(), key=lambda x: -x[1])

        if all([snp['gene'] == 'null' for snp in self.analyzed_variants]):
            gene_snp_rate = []

        fname = os.path.join(self.output_dir, '{}_top_genes_var_rate.svg'.format(self.study_name))
        genes, rates = zip(*gene_snp_rate[:100])  # top 100
        self.generate_bar_chart([genes, rates], 'gene_var_rate', fname)

        if write_summary_pkl:
            self.summary = {
                'study_name': self.study_name,
                'Counts summary': {
                    'Total number of variants': len(snp_set),
                    'Total number of affected genes': len(gene_snp_rate),
                    'Total number of affected proteins': len(set(snp['ref_prot_acc'] for snp in self.analyzed_variants if snp['ref_prot_acc'] != 'null')),
                },
                'Chromosome Level': snps_per_chroms,
                'Gene Level': gene_snp_rate,
            }
            fname = os.path.join(self.output_dir, 'summary.pkl')
            pickle.dump(self.summary, open(fname, 'wb'))

        else:
            return snps_per_chroms, gene_snp_rate

    def compute_summary_stats(self):
        # overall summary
        snp_set = set([snp['hg38_name'] for snp in self.analyzed_variants])
        prot_counts = {
            'alt': len(set(snp['alt_prot_acc'] for snp in self.analyzed_variants if 'IP_' in snp['alt_prot_acc'])),
            'iso': len(set(snp['alt_prot_acc'] for snp in self.analyzed_variants if 'II_' in snp['alt_prot_acc'])),
            'ref': len(set(snp['ref_prot_acc'] for snp in self.analyzed_variants if snp['ref_prot_acc'] != 'null'))
        }

        # chrom level & gene level
        snps_per_chroms, gene_snp_rate = self.compute_chrom_gene_level_stats()

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
        gene_snps_grp = sorted(self.analyzed_variants, key=lambda x: x['gene'])
        if len(self.analyzed_variants) == len([x for x in self.analyzed_variants if x['gene'] == 'null']):
            alt_snps_stats = {'All nulls': 'No gene consequences for the submitted variants.'}
        else:
            alt_grp_snp = sorted(self.analyzed_variants, key = lambda x: x['alt_prot_acc'])
            highest_impact_ann = {'alt_highest': [snp['alt_prot_acc'] for snp in self.analyzed_variants if snp['alt_max_impact'] > snp['ref_max_impact']], 'ref_highest': [snp['ref_prot_acc'] for snp in self.analyzed_variants if snp['ref_max_impact'] > snp['alt_max_impact']]}
            gene_snps = {gene: len(list(grp)) for gene, grp in itt.groupby(gene_snps_grp, key = lambda x: x['gene'])}
            for snp in alt_grp_snp:
                snp.update( {'total_gene_snp': gene_snps[snp['gene']]} )
            alt_snps_stats = {gene_alt: self.count_altsnp_stats(list(grp)) for gene_alt, grp in itt.groupby(alt_grp_snp, key = lambda x: x['gene'] + ' | ' + x['alt_prot_acc']) if (gene_alt.split(' | ')[1] in highest_impact_ann['alt_highest']) and (gene_alt.split(' | ')[1] != 'null')}


        self.summary = {
	        'study_name': self.study_name,
            'Counts summary': {
                'Total number of variants': len(snp_set),
                'Total number of affected genes': len(gene_snps),
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
                'Impact Annotation': impact_ann,
                'Fold Change': fc,
            },
            'Mutational hotspots on altORFs': alt_snps_stats,
        }

        if len(self.analyzed_variants) != len([x for x in self.analyzed_variants if x['gene'] == 'null']):
            fname = os.path.join(self.output_dir, '{}_hotspots_barchart.svg'.format(self.study_name))
            self.generate_bar_chart(alt_snps_stats, 'hotspots_bar', fname)

        fname = os.path.join(self.output_dir, 'summary.pkl')
        pickle.dump(self.summary, open(fname, 'wb'))

    def generate_bar_chart(self, data, chart_type, fname):
        if chart_type == 'gene_var_rate':
            genes, rates = data
            fig, ax = plt.subplots()
            ax.bar(range(1, len(genes) + 1), rates, color='#bde2f8')
            plt.xticks(range(1, len(genes) + 1), genes, rotation='vertical')
            plt.title('Mutations per Kb per gene')
            plt.xlabel('Genes')
            plt.ylabel('SNPs per Kb')
            plt.savefig(fname)
            plt.show()

        if chart_type == 'fold_change':
            data = [data[i] for i in range(1, 4)]
            fig, ax = plt.subplots()
            ax.bar(range(1, 4), data, color=['#e6f2ff', '#66b3ff', '#004d99'])
            ax.set_facecolor('dimgrey')
            plt.xticks(range(1, 4), ['Low', 'Medium', 'High'])
            plt.title('Fold-change gained from deeper annotation per impact')
            plt.xlabel('Impact Levels')
            plt.ylabel('Fold-change gained from deeper genome annotation')
            plt.savefig(fname)
            plt.show()

        if chart_type == 'stacked_impact':
            ref_impacts, alt_impacts = data
            fig, ax = plt.subplots()
            ax.bar(range(1, 4), ref_impacts, label='Canonical ORFs', color='#39ac73')
            ax.bar(range(1, 4), alt_impacts, bottom=ref_impacts, label='Alternative ORFs', color='#9fdfbf')
            ax.set_facecolor('dimgrey')
            plt.xticks(range(1, 4), ['Low', 'Medium', 'High'])
            plt.xlabel('Impact Levels')
            plt.ylabel('Count of SNPs')
            plt.legend()
            plt.savefig(fname)
            plt.show()

        if chart_type == 'hotspots_bar':
            bins = [(0. + (n - 1) * (1. / 30), 0. + n * (1. / 30)) for n in list(range(1, 31))]
            bin_labels = [' - '.join(['{:.2f}'.format(round(x, 2)) for x in left_right]) for left_right in bins]
            altorfs_per_bin = {cat: {n: [] for n in bin_labels} for cat in ['one_snp', 'one_ten', 'over_ten']}
            for gene_alt in data:
                for left, right, label in zip([x[0] for x in bins], [x[1] for x in bins], bin_labels):
                    freq = data[gene_alt]['total_portion_gene']
                    snps = data[gene_alt]['cnt_snps']
                    mean_impact = data[gene_alt]['mean_impacts']
                    if (freq > left) and (freq <= right):
                        if snps == 1:
                            altorfs_per_bin['one_snp'][label].append(mean_impact)
                        elif snps < 10:
                            altorfs_per_bin['one_ten'][label].append(mean_impact)
                        else:
                            altorfs_per_bin['over_ten'][label].append(mean_impact)
            altorf_counts = {cat: {n: len(altorfs_per_bin[cat][n]) for n in bin_labels} for cat in ['one_snp', 'one_ten', 'over_ten']}
            colors = ['#f3a6fc', '#cf5fe3', '#7b198c']
            fig, ax = plt.subplots()
            ax.bar(list(range(1, 31, 1)), list(altorf_counts['one_snp'].values()), color = colors[0], label = 'One SNP')
            ax.bar(list(range(1, 31, 1)), list(altorf_counts['one_ten'].values()), bottom = list(altorf_counts['one_snp'].values()), color = colors[1], label = 'Less than 10 SNPs')
            ax.bar(list(range(1, 31, 1)), list(altorf_counts['over_ten'].values()), bottom = [altorf_counts['one_snp'][n] + altorf_counts['one_ten'][n] for n in bin_labels], color = colors[2], label = '10 or more SNPs')
            ax.set_facecolor('dimgrey')
            plt.xticks(list(range(1, 31)), bin_labels, rotation = 90)
            plt.title('Mutational hotspots on altORFs')
            plt.xlabel("Portion of a gene's SNPs within altORF")
            plt.ylabel('Count of altORFs')
            plt.legend()
            plt.savefig(fname)
            plt.show()

    def count_altsnp_stats(self, snps_dict):
        cnt_snps = len(set(snp['hg38_name'] for snp in snps_dict))
        impacts = list(snp['alt_max_impact'] for snp in snps_dict)
        hi_impacts = list(snp['alt_max_impact'] for snp in snps_dict if snp['alt_max_impact'] >= snp['ref_max_impact'])
        tot_gene_snps = list(snp['total_gene_snp'] for snp in snps_dict)[0]
        return {'cnt_snps': cnt_snps, 'impacts': impacts, 'mean_impacts': np.mean([snp['alt_max_impact'] for snp in snps_dict]), 'highest_impacts': hi_impacts,
                'total_portion_gene': cnt_snps / tot_gene_snps, 'hi_portion_gene': len(hi_impacts) / tot_gene_snps,
                'hi_portion_alt': len(hi_impacts) / len(impacts), 'score': (cnt_snps ** np.mean([snp['alt_max_impact'] for snp in snps_dict]))}

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
            print('All genes are null!')
            #raise Exception('All genes are null!')

        print('All variants analyzed')
        return analyzed_variants

    def list_annOnePerLine_files(self):
        annOnePerLine_files = []
        _dir = self.opv.vcf.results_dir
        if 'OP_' in self.opv.annotation:
            _dir = self.opv.vcf.vcf_splits_dir

        for f in os.listdir(_dir):
            fpath = os.path.join(_dir, f)
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
                var_name = '_'.join([line['CHROM'], line['POS'], line['REF'], line['ALT'], line['ANN[*].GENE']])
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
