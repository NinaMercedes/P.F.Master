#!/usr/bin/env python

import sys
import argparse
import os
import fastq2matrix as fm

def get_step_num(prefix_path):
    files = {
        f"{prefix_path}.mkdup.bam.bai": 1,
        f"{prefix_path}.bqsr.bam.bai": 2,
        f"{prefix_path}.bqsr.cram.crai": 2,
        f"{prefix_path}.mkdup.cram.crai": 2,
        f"{prefix_path}.g.vcf.gz.validated": 3
    }
    step = 0
    for f in files:
        if os.path.isfile(f):
            step = files[f]
            sys.stderr.write(f"Found {f}\n")
    return step

def main_trim(args):
    args.prefix_path = os.path.join(args.tmp_dir, args.prefix)
    if args.single:
        trimmed_file = f"{args.prefix_path}_trimmed.fq"
        fm.run_cmd(f"trimmomatic SE -phred33 {args.read1} {trimmed_file} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
    else:
        fm.run_cmd(f"trimmomatic PE -phred33 {args.read1} {args.read2} -baseout {args.prefix_path} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
    args.trimmed = True

def main_map(args):
    args.prefix_path = os.path.join(args.tmp_dir, args.prefix)
    args.step = get_step_num(args.prefix_path)

    # Determine which reads to map
    if getattr(args, "trimmed", False):
        if args.single:
            reads = f"{args.prefix_path}_trimmed.fq"
        else:
            reads = f"{args.prefix_path}_1P {args.prefix_path}_2P"
    else:
        if args.single:
            reads = args.read1
        else:
            reads = f"{args.read1} {args.read2}"

    drop_unmapped = "-F 4" if getattr(args, "drop_unmapped", False) else ""

    # Mapping command
    fm.run_cmd(f"{args.mapper} mem -t {args.threads} -R \"@RG\\tID:{args.prefix}\\tSM:{args.prefix}\\tPL:Illumina\" {args.ref} {reads} | "
               f"samtools view -@ {args.threads} {drop_unmapped} -b - | "
               f"samtools fixmate -@ {args.threads} -m - - | "
               f"samtools sort -@ {args.threads} - | "
               f"samtools markdup -@ {args.threads} - {args.prefix_path}.mkdup.bam -")

    # Remove trimmed files
    if getattr(args, "trimmed", False):
        if args.single:
            fm.run_cmd(f"rm -f {reads}")
        else:
            fm.run_cmd(f"rm -f {args.prefix_path}_1P {args.prefix_path}_2P {args.prefix_path}_1U {args.prefix_path}_2U")

    # Index and flagstat
    fm.run_cmd(f"samtools index -@ {args.threads} {args.prefix_path}.mkdup.bam")
    fm.run_cmd(f"samtools flagstat -@ {args.threads} {args.prefix_path}.mkdup.bam > {args.prefix_path}.mkdup.bamstats")

    # BQSR
    if args.bqsr_vcf and (args.redo or args.step < 2):
        for vcf in args.bqsr_vcf.split(","):
            fm.tabix_vcf(vcf)
        bqsr_sites = " ".join([f"--known-sites {s}" for s in args.bqsr_vcf.split(",")])
        fm.run_cmd(f"gatk BaseRecalibrator -R {args.ref} -I {args.prefix_path}.mkdup.bam {bqsr_sites} -O {args.prefix_path}.recal_data.table")
        fm.run_cmd(f"gatk ApplyBQSR -R {args.ref} -I {args.prefix_path}.mkdup.bam --bqsr-recal-file {args.prefix_path}.recal_data.table -O {args.prefix_path}.bqsr.bam")
        fm.run_cmd(f"samtools index -@ {args.threads} {args.prefix_path}.bqsr.bam")
        fm.run_cmd(f"samtools flagstat -@ {args.threads} {args.prefix_path}.bqsr.bam > {args.prefix_path}.bqsr.bamstats")
        fm.run_cmd(f"rm {args.prefix_path}.mkdup.bam*")
        args.bam = f"{args.prefix_path}.bqsr.bam"
    else:
        args.bam = f"{args.prefix_path}.mkdup.bam"

def main_gatk(args):
    if not args.prefix:
        args.prefix = args.bam.replace(".bam", "")
    fm.run_cmd(f"gatk HaplotypeCaller -I {args.bam} -R {args.ref} -O {args.prefix}.g.vcf.gz -ERC {args.erc}")
    fm.run_cmd(f"gatk ValidateVariants -V {args.prefix}.g.vcf.gz -gvcf -R {args.ref} && touch {args.prefix}.g.vcf.gz.validated")

def convert_to_cram(bam_file, ref_file, threads):
    cram_file = bam_file.replace(".bam", ".cram")
    fm.run_cmd(f"samtools view -@ {threads} -C {bam_file} -o {cram_file} -T {ref_file}")
    fm.run_cmd(f"samtools index {cram_file}")
    fm.run_cmd(f"rm {bam_file} {bam_file}.bai")
    return cram_file

def main_all(args):
    args.prefix_path = os.path.join(args.tmp_dir, args.prefix)
    fm.create_seq_dict(args.ref)
    if args.mapper == "bwa":
        fm.bwa_index(args.ref)
    elif args.mapper == "bwa-mem2":
        fm.bwa2_index(args.ref)
    fm.faidx(args.ref)

    args.step = get_step_num(args.prefix_path)

    if not args.single and not args.read2:
        sys.stderr.write("Second read is missing for paired-end run. Exiting.\n")
        sys.exit(1)

    sys.stderr.write(f"Starting at step {args.step+1}\n")
    if args.trim and (args.redo or args.step < 1):
        sys.stderr.write("Running trimming step...\n")
        main_trim(args)
    else:
        sys.stderr.write("Skipping trimming step (use --trim to enable)\n")
        args.trimmed = False

    main_map(args)

    if args.redo or args.step < 3:
        sys.stderr.write(f"Using {args.bam} as the bam file\n")
        main_gatk(args)

    if args.cram:
        if args.redo or args.step < 4:
            args.bam = convert_to_cram(args.bam, args.ref, args.threads)

# ---------------- Argument parsing ----------------
parser = argparse.ArgumentParser(description='fastq2matrix pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

# All
parser_sub = subparsers.add_parser('all', help='Run full pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1','-1', help='First read file', required=True)
parser_sub.add_argument('--read2','-2', help='Second read file')
parser_sub.add_argument('--prefix','-p', help='Sample prefix', required=True)
parser_sub.add_argument('--ref','-r', help='Reference genome', required=True)
parser_sub.add_argument('--threads','-t', default=4, type=int, help='Number of threads')
parser_sub.add_argument('--bqsr-vcf','-q', help='VCF file used for BQSR')
parser_sub.add_argument('--erc', default="GVCF", choices=["GVCF","BP_RESOLUTION"], help='ERC type for GATK')
parser_sub.add_argument('--redo', action="store_true", help='Redo all steps')
parser_sub.add_argument('--single', action="store_true", help='Single-end reads')
parser_sub.add_argument('--cram', action="store_true", help='Convert BAM to CRAM')
parser_sub.add_argument('--trim', action='store_true', help='Run trimming before mapping')
parser_sub.add_argument('--mapper','-m', default="bwa", choices=["bwa","bwa-mem2"], help='Mapper to use')
parser_sub.add_argument('--tmp-dir', default=".", help='Temporary working directory')
parser_sub.set_defaults(func=main_all)

# Trim
parser_sub = subparsers.add_parser('trim', help='Trim reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--prefix','-p', help='Sample prefix', required=True)
parser_sub.set_defaults(func=main_trim)

# Map
parser_sub = subparsers.add_parser('map', help='Map reads', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1','-1', help='First read file', required=True)
parser_sub.add_argument('--read2','-2', help='Second read file', required=True)
parser_sub.add_argument('--prefix','-p', help='Sample prefix', required=True)
parser_sub.add_argument('--ref','-r', help='Reference genome', required=True)
parser_sub.add_argument('--threads','-t', default=4, type=int, help='Number of threads')
parser_sub.add_argument('--bqsr-vcf','-q', help='VCF file used for BQSR')
parser_sub.add_argument('--redo', action="store_true", help='Redo all steps')
parser_sub.add_argument('--single', action="store_true", help='Single-end reads')
parser_sub.add_argument('--cram', action="store_true", help='Convert BAM to CRAM')
parser_sub.add_argument('--mapper','-m', default="bwa", choices=["bwa","bwa-mem2"], help='Mapper to use')
parser_sub.set_defaults(func=main_map)

# GATK
parser_sub = subparsers.add_parser('gatk', help='Run GATK variant calling', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--bam','-b', help='Input BAM file', required=True)
parser_sub.add_argument('--ref','-r', help='Reference genome', required=True)
parser_sub.add_argument('--prefix','-p', help='Sample prefix')
parser_sub.add_argument('--erc', default="GVCF", choices=["GVCF","BP_RESOLUTION"], help='ERC type for GATK')
parser_sub.add_argument('--threads','-t', default=4, type=int, help='Number of threads')
parser_sub.set_defaults(func=main_gatk)

# ---------------- Parse args ----------------
args = parser.parse_args()
if vars(args) == {}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
