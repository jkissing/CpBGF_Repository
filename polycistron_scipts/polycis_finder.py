#!/usr/bin/env python3
import sys
import os
import subprocess
import argparse
import tempfile
from datetime import datetime
from collections import defaultdict
def check_and_install_modules(modules):
    import importlib
    missing_modules = []
    for module in modules:
        try:
            importlib.import_module(module)
        except ImportError:
            missing_modules.append(module)
    if missing_modules:
        print(f"The following Python modules are missing: {', '.join(missing_modules)}")
        install = input("Do you want to install them now? [y/n]: ").strip().lower()
        if install == 'y':
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install"] + missing_modules)
            except subprocess.CalledProcessError:
                print("Failed to install required Python modules. Please install them manually.", file=sys.stderr)
                sys.exit(1)
        else:
            print("Cannot proceed without installing the required Python modules.", file=sys.stderr)
            sys.exit(1)
required_modules = ['pybedtools', 'tqdm']
check_and_install_modules(required_modules)
import pybedtools
from tqdm import tqdm
def check_bedtools():
    """Check if bedtools is installed and accessible."""
    try:
        subprocess.run(['bedtools', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("Error: bedtools is not installed or not found in PATH.", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: bedtools is not installed or not found in PATH.", file=sys.stderr)
        sys.exit(1)

def sort_bedfile(input_bed):
    """Sort a BED file using pybedtools."""
    try:
        bed = pybedtools.BedTool(input_bed)
        sorted_bed = bed.sort()
        return sorted_bed
    except Exception as e:
        print(f"Error sorting BED file {input_bed}: {e}", file=sys.stderr)
        sys.exit(1)

def parse_cds(sorted_cds_bed):
    """
    Parse the sorted CDS BED file.
    Returns:
        gene_to_cds: dict mapping gene_name to set of CDS identifiers
    """
    gene_to_cds = defaultdict(set)
    try:
        with tqdm(sorted_cds_bed, desc="Parsing CDS annotations", unit=" lines") as iter_bed:
            for interval in iter_bed:
                if not interval.fields:
                    continue
                if len(interval.fields) < 6:
                    print(f"Warning: Skipping malformed CDS line: {str(interval)}", file=sys.stderr)
                    continue
                cds_info = interval.name
                try:
                    cds_name, gene_name = cds_info.split(';')
                except ValueError:
                    print(f"Warning: Unable to split CDS and gene name in line: {str(interval)}", file=sys.stderr)
                    continue
                gene_to_cds[gene_name].add(cds_name)
    except Exception as e:
        print(f"Error parsing CDS BED file: {e}", file=sys.stderr)
        sys.exit(1)
    return gene_to_cds

def parse_transcript_bed(sorted_transcript_bed):
    """
    Parse the sorted transcript BED file.
    Returns:
        transcript_info: dict mapping transcript_name to its (chrom, start, end, strand)
    """
    transcript_info = {}
    try:
        with tqdm(sorted_transcript_bed, desc="Parsing transcript models", unit=" lines") as iter_bed:
            for interval in iter_bed:
                if not interval.fields:
                    continue
                if len(interval.fields) < 6:
                    print(f"Warning: Skipping malformed transcript line: {str(interval)}", file=sys.stderr)
                    continue
                chrom = interval.chrom
                start = interval.start
                end = interval.end
                transcript_name = interval.name
                strand = interval.strand
                transcript_info[transcript_name] = (chrom, start, end, strand)
    except Exception as e:
        print(f"Error parsing transcript BED file: {e}", file=sys.stderr)
        sys.exit(1)
    return transcript_info

def run_bedtools_intersect(sorted_transcript, sorted_cds):
    """Run bedtools intersect to find CDS fully covered by transcripts on the same strand."""
    try:
        intersect = sorted_transcript.intersect(sorted_cds, f=1.0, s=True, wa=True, wb=True)
        return intersect
    except Exception as e:
        print(f"Error running bedtools intersect: {e}", file=sys.stderr)
        sys.exit(1)

def process_intersection(intersect, gene_to_cds):
    """
    Process the bedtools intersect output.
    Returns:
        transcript_to_genes: dict mapping transcript_name to set of genes fully covered
    """
    transcript_to_covered_cds = defaultdict(lambda: defaultdict(set))
    try:
        with tqdm(intersect, desc="Processing intersections", unit=" lines") as iter_intersect:
            for interval in iter_intersect:
                if not interval.fields:
                    continue
                fields = interval.fields
                if len(fields) < 12:
                    print(f"Warning: Skipping malformed intersect line: {str(interval)}", file=sys.stderr)
                    continue
                t_name = fields[3]
                b_info = fields[9]
                try:
                    cds_name, gene_name = b_info.split(';')
                except ValueError:
                    print(f"Warning: Unable to split CDS and gene name in intersect line: {str(interval)}", file=sys.stderr)
                    continue
                transcript_to_covered_cds[t_name][gene_name].add(cds_name)
    except Exception as e:
        print(f"Error processing intersection data: {e}", file=sys.stderr)
        sys.exit(1)
    
    transcript_to_genes = defaultdict(set)
    for transcript, genes in transcript_to_covered_cds.items():
        for gene, covered_cds in genes.items():
            total_cds = gene_to_cds.get(gene, set())
            if total_cds and covered_cds == total_cds:
                transcript_to_genes[transcript].add(gene)
    return transcript_to_genes

def generate_output_filename(input1, input2):
    """Generate the output filename based on input filenames and current date."""
    base1 = os.path.splitext(os.path.basename(input1))[0]
    base2 = os.path.splitext(os.path.basename(input2))[0]
    date_str = datetime.now().strftime("%Y%m%d")
    output_filename = f"{base1}_{base2}_output_{date_str}.tsv"
    return output_filename

def write_output(output_filename, transcript_info, transcript_to_genes):
    """Write the output TSV file."""
    try:
        with open(output_filename, 'w') as out:
            header = ["transcript_model", "chr", "start", "end", "strand", "genes_covered", "gene_names"]
            out.write('\t'.join(header) + '\n')
            for transcript, info in tqdm(transcript_info.items(), desc="Writing output", unit=" transcripts"):
                chrom, start, end, strand = info
                genes = transcript_to_genes.get(transcript, set())
                genes_covered = len(genes)
                gene_names = ';'.join(sorted(genes)) if genes else ""
                row = [transcript, chrom, str(start), str(end), strand, str(genes_covered), gene_names]
                out.write('\t'.join(row) + '\n')
    except Exception as e:
        print(f"Error writing output file {output_filename}: {e}", file=sys.stderr)
        sys.exit(1)

def generate_output_filename(input1, input2):
    """Generate the output filename based on input filenames and current date."""
    base1 = os.path.splitext(os.path.basename(input1))[0]
    base2 = os.path.splitext(os.path.basename(input2))[0]
    date_str = datetime.now().strftime("%Y%m%d")
    output_filename = f"{base1}_{base2}_output_{date_str}.tsv"
    return output_filename

def main():
    parser = argparse.ArgumentParser(
        description="Generate a report on how many genes' CDS are fully covered by each transcript model.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("transcript_bed", help="Input transcript BED file (bed format).")
    parser.add_argument("cds_bed", help="Input CDS annotation BED file (bed format).")
    args = parser.parse_args()

    check_bedtools()

    with tempfile.TemporaryDirectory() as tmpdir:
        print("Sorting transcript BED file...")
        sorted_transcript = sort_bedfile(args.transcript_bed)
        print("Sorting CDS BED file...")
        sorted_cds = sort_bedfile(args.cds_bed)

        print("Parsing CDS annotations...")
        gene_to_cds = parse_cds(sorted_cds)

        print("Running bedtools intersect...")
        intersect = run_bedtools_intersect(sorted_transcript, sorted_cds)

        print("Parsing transcript models...")
        transcript_info = parse_transcript_bed(sorted_transcript)

        print("Processing intersections to determine fully covered genes...")
        transcript_to_genes = process_intersection(intersect, gene_to_cds)

        output_filename = generate_output_filename(args.transcript_bed, args.cds_bed)
        print(f"Writing output to {output_filename}...")
        write_output(output_filename, transcript_info, transcript_to_genes)

        print("Processing complete.")

if __name__ == "__main__":
    main()