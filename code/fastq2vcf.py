#! /usr/bin/env python

import sys
import argparse
import fastq2matrix as fm
import os.path


def get_step_num(prefix):
    files = {
        f"{prefix}.mkdup.bam.bai": 1,
        f"{prefix}.bqsr.bam.bai": 2,
        f"{prefix}.bqsr.cram.crai": 2,
        f"{prefix}.mkdup.cram.crai": 2,
        f"{prefix}.g.vcf.gz.validated": 3
    }
    step = 0
    for f in files:
        if os.path.isfile(f):
            step = files[f]
            sys.stderr.write(f"Found {f}\n")
    return step


def main_trim(args):
    args.prefix_path = args.tmp_dir + "/" + args.prefix
    if args.single:
        fm.run_cmd(
            "trimmomatic SE -phred33 -threads %(threads)s %(read1)s %(tmp_dir)s/%(prefix)s_trimmed.fq "
            "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(args)
        )
    else:
        fm.run_cmd(
            "trimmomatic PE -phred33 -threads %(threads)s %(read1)s %(read2)s "
            "-baseout %(tmp_dir)s/%(prefix)s LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36" % vars(args)
        )


def main_map(args):
    if args.mapper == "bwa":
        fm.bwa_index(args.ref)
    elif args.mapper == "bwa-mem2":
        fm.bwa2_index(args.ref)

    args.prefix_path = args.tmp_dir + "/" + args.prefix
    args.step = get_step_num(args.prefix_path)

    # Determine which reads to use
    if args.trimmed and args.single:
        args.reads = "%(tmp_dir)s/%(prefix)s_trimmed.fq" % vars(args)
    elif args.trimmed and not args.single:
        args.reads = "%(tmp_dir)s/%(prefix)s_1P %(tmp_dir)s/%(prefix)s_2P" % vars(args)
    elif not args.trimmed and not args.single:
        args.reads = "%(read1)s %(read2)s" % vars(args)
    elif not args.trimmed and args.single:
        args.reads = "%(read1)s" % vars(args)

    if args.redo or args.step < 1:
        args.drop_unmapped = "-F 4 " if args.drop_unmapped else ""
        fm.run_cmd(
            "%(mapper)s mem -t %(threads)s -R \"@RG\\tID:%(prefix)s\\tSM:%(prefix)s\\tPL:Illumina\" "
            "%(ref)s %(reads)s | samtools view -@ %(threads)s %(drop_unmapped)s -b - | "
            "samtools fixmate -@ %(threads)s -m - - | samtools sort -@ %(threads)s - | "
            "samtools markdup -@ %(threads)s - %(tmp_dir)s/%(prefix)s.mkdup.bam -" % vars(args)
        )

        if args.trimmed and args.single:
            fm.run_cmd("rm %(reads)s" % vars(args))
        if args.trimmed and not args.single:
            fm.run_cmd(
                "rm %(tmp_dir)s/%(prefix)s_1P %(tmp_dir)s/%(prefix)s_2P "
                "%(tmp_dir)s/%(prefix)s_1U %(tmp_dir)s/%(prefix)s_2U" % vars(args)
            )

        fm.run_cmd("samtools index -@ %(threads)s %(tmp_dir)s/%(prefix)s.mkdup.bam" % vars(args))
        args.bam = "%(tmp_dir)s/%(prefix)s.mkdup.bam" % vars(args)

    if args.bqsr_vcf and (args.redo or args.step < 2):
        for vcf in args.bqsr_vcf.split(","):
            fm.tabix_vcf(vcf)
        args.bqsr_vcf = " ".join(["--known-sites %s" % s for s in args.bqsr_vcf.split(",")])
        fm.run_cmd(
            "gatk BaseRecalibrator -R %(ref)s -I %(tmp_dir)s/%(prefix)s.mkdup.bam "
            "%(bqsr_vcf)s -O %(tmp_dir)s/%(prefix)s.recal_data.table" % vars(args)
        )
        fm.run_cmd(
            "gatk ApplyBQSR -R %(ref)s -I %(tmp_dir)s/%(prefix)s.mkdup.bam "
            "--bqsr-recal-file %(tmp_dir)s/%(prefix)s.recal_data.table -O %(tmp_dir)s/%(prefix)s.bqsr.bam" % vars(args)
        )
        fm.run_cmd("samtools index -@ %(threads)s %(tmp_dir)s/%(prefix)s.bqsr.bam" % vars(args))
        fm.run_cmd("rm %(tmp_dir)s/%(prefix)s.mkdup.bam*" % vars(args))
        args.bam = "%(tmp_dir)s/%(prefix)s.bqsr.bam" % vars(args)

    if args.cram:
        cram_file = args.bam.replace(".bam", ".cram")
        if args.redo or not os.path.isfile(cram_file):
            sys.stderr.write("Converting to CRAM\n")
            convert_to_cram(args.bam, args.ref, args.threads)
        args.bam = cram_file

    if args.bam_qc:
        if args.redo or not os.path.isfile(args.bam + ".flagstat.txt") or not os.path.isfile(args.bam + ".genomecov.txt"):
            sys.stderr.write("Creating coverage and flagstat files\n")
            bam_qc(args)


def main_gatk(args):
    args.prefix_path = args.tmp_dir + "/" + args.prefix
    if not args.prefix:
        args.prefix = args.bam.replace(".bam", "")
    args.regions = "-L %s" % args.gatk_bed if args.gatk_bed else ""
    fm.run_cmd(
        "gatk HaplotypeCaller -I %(bam)s -R %(ref)s -O %(tmp_dir)s/%(prefix)s.g.vcf.gz "
        "-ERC %(erc)s %(hc_options)s %(regions)s" % vars(args)
    )
    fm.run_cmd(
        "gatk ValidateVariants -V %(tmp_dir)s/%(prefix)s.g.vcf.gz -gvcf -R %(ref)s %(regions)s "
        "&& touch %(tmp_dir)s/%(prefix)s.g.vcf.gz.validated" % vars(args)
    )


def convert_to_cram(bam_file, ref_file, threads):
    cram_file = bam_file.replace(".bam", ".cram")
    fm.run_cmd("samtools view -@ %s -C %s -o %s -T %s" % (threads, bam_file, cram_file, ref_file))
    fm.run_cmd("samtools index %s" % cram_file)
    fm.run_cmd("rm %s %s.bai" % (bam_file, bam_file))
    return cram_file


def bam_qc(args):
    fm.run_cmd("samtools flagstat -@ %(threads)s %(bam)s > %(bam)s.flagstat.txt" % vars(args))
    fm.run_cmd("samtools stats -@ %(threads)s %(bam)s > %(bam)s.stats.txt" % vars(args))
    fm.run_cmd("bedtools genomecov -ibam %(bam)s > %(bam)s.genomecov.txt" % vars(args))


def main_all(args):
    fm.create_seq_dict(args.ref)
    fm.faidx(args.ref)

    args.prefix_path = args.tmp_dir + "/" + args.prefix
    args.step = get_step_num(args.prefix_path)

    if not args.read2 and not args.single:
        sys.stderr.write("Second read is not provided, please check... Exiting!\n")
        quit()

    args.bam = (
        args.tmp_dir + "/" + args.prefix + ".bqsr.bam"
        if args.bqsr_vcf
        else args.tmp_dir + "/" + args.prefix + ".mkdup.bam"
    )

    sys.stderr.write("Starting at step %s\n" % (args.step + 1))

    # Optional trimming step
    if args.trim and (args.redo or args.step < 1):
        sys.stderr.write("Running trimming step\n")
        main_trim(args)
        args.trimmed = True
    else:
        args.trimmed = False

    # Mapping step
    if args.redo or args.step < 2:
        main_map(args)

    # Variant calling
    if (args.redo or args.step < 3) and not args.no_variant_calling:
        sys.stderr.write("Using %(bam)s as the bam file\n" % vars(args))
        main_gatk(args)

    # Move results
    if os.path.abspath(args.storage_dir) != os.getcwd():
        fm.run_cmd("mv %(prefix_path)s* %(storage_dir)s/" % vars(args))


parser = argparse.ArgumentParser(description='fastq2matrix pipeline', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

# --- ALL ---
parser_sub = subparsers.add_parser('all', help='Run full pipeline: optional trim → map → GATK',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1', '-1', help='First read file', required=True)
parser_sub.add_argument('--read2', '-2', help='Second read file')
parser_sub.add_argument('--prefix', '-p', help='Sample prefix for all results generated', required=True)
parser_sub.add_argument('--ref', '-r', help='Reference genome FASTA', required=True)
parser_sub.add_argument('--threads', '-t', default=4, help='Number of threads')
parser_sub.add_argument('--mapper', '-m', default="bwa", choices=["bwa", "bwa-mem2"], type=str, help='Mapping program to use')
parser_sub.add_argument('--bqsr-vcf', '-q', help='VCF file(s) used for BQSR')
parser_sub.add_argument('--erc', default="GVCF", choices=["GVCF", "BP_RESOLUTION"], help='Choose ERC type for GATK')
parser_sub.add_argument('--hc-options', default="", type=str, help='Additional HaplotypeCaller options')
parser_sub.add_argument('--redo', action="store_true", help='Redo all steps')
parser_sub.add_argument('--single', action="store_true", help='Single-end reads')
parser_sub.add_argument('--cram', action="store_true", help='Convert BAM to CRAM after mapping')
parser_sub.add_argument('--bam-qc', action="store_true", help='Generate flagstat and coverage reports')
parser_sub.add_argument('--tmp-dir', default=".", type=str, help='Temporary working directory')
parser_sub.add_argument('--storage-dir', default=".", type=str, help='Final results directory')
parser_sub.add_argument('--gatk-bed', type=str, help='Target BED regions for GATK')
parser_sub.add_argument('--drop-unmapped', action="store_true", help='Drop unmapped reads during mapping')
parser_sub.add_argument('--no-variant-calling', action="store_true", help="Skip variant calling")
parser_sub.add_argument('--trim', action="store_true", help='Perform trimming with Trimmomatic before mapping')
parser_sub.set_defaults(func=main_all)

# --- TRIM ---
parser_sub = subparsers.add_parser('trim', help='Trim reads using Trimmomatic',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--prefix', '-p', required=True)
parser_sub.add_argument('--read1', '-1', required=True)
parser_sub.add_argument('--read2', '-2')
parser_sub.add_argument('--single', action="store_true")
parser_sub.add_argument('--threads', '-t', default=4)
parser_sub.add_argument('--tmp-dir', default=".", type=str)
parser_sub.set_defaults(func=main_trim)

# --- MAP ---
parser_sub = subparsers.add_parser('map', help='Map reads using BWA/BWA-MEM2',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--read1', '-1', required=True)
parser_sub.add_argument('--read2', '-2')
parser_sub.add_argument('--prefix', '-p', required=True)
parser_sub.add_argument('--ref', '-r', required=True)
parser_sub.add_argument('--threads', '-t', default=4)
parser_sub.add_argument('--bqsr-vcf', '-q')
parser_sub.add_argument('--mapper', '-m', default="bwa", choices=["bwa", "bwa-mem2"])
parser_sub.add_argument('--redo', action="store_true")
parser_sub.add_argument('--tmp-dir', default=".")
parser_sub.add_argument('--storage-dir', default=".")
parser_sub.add_argument('--drop-unmapped', action="store_true")
parser_sub.add_argument('--single', action="store_true")
parser_sub.add_argument('--cram', action="store_true")
parser_sub.add_argument('--bam-qc', action="store_true")
parser_sub.set_defaults(func=main_map)

# --- GATK ---
parser_sub = subparsers.add_parser('gatk', help='Run GATK HaplotypeCaller',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--bam', '-b', required=True)
parser_sub.add_argument('--ref', '-r', required=True)
parser_sub.add_argument('--prefix', '-p')
parser_sub.add_argument('--erc', default="GVCF", choices=["GVCF", "BP_RESOLUTION"])
parser_sub.add_argument('--hc-options', default="", type=str)
parser_sub.add_argument('--threads', '-t', default=4)
parser_sub.add_argument('--tmp-dir', default=".")
parser_sub.add_argument('--storage-dir', default=".")
parser_sub.add_argument('--gatk-bed', type=str)
parser_sub.set_defaults(func=main_gatk)


args = parser.parse_args()
if vars(args) == {}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
