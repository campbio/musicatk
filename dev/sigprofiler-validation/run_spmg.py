"""Count the musicatk example variants with SigProfilerMatrixGenerator.

One-off validation for issue #218; not part of the package. See README.md
in this folder for how to set up the environment and rerun.

Usage, from the package root:
    python dev/sigprofiler-validation/run_spmg.py inst/extdata <out_dir>
    python dev/sigprofiler-validation/run_spmg.py --dir <input_dir> \
        <project> <genome> <out_dir>

The second form runs SPMG on one folder of VCF or MAF files that is already
prepared (compare_tcga.R and compare_pcawg.R use it). The first form runs
the checks on the bundled example files:

The four VCFs are copied with only their PASS records, because SPMG ignores
the FILTER column while musicatk keeps PASS VCF records by default. The LUSC
MAF is copied unchanged; neither tool filters MAF records. SPMG runs on
GRCh38 for four projects:

- "vcf": the four VCFs
- "maf": the LUSC MAF
- "insfix": the four MAF insertions written as VCF records anchored on the
  base at Start_Position (REF is that hg38 base), to check SPMG's insertion
  bins without its MAF coordinate shift
- "synthetic": a run of three adjacent SBS and a complex indel, which the
  example files don't have
"""

import os
import shutil
import sys

from SigProfilerMatrixGenerator.scripts import (
    SigProfilerMatrixGeneratorFunc as matGen,
)

VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

# sample: [(chrom, pos, ref, alt)]; REF bases are hg38
EXTRA = {
    "insfix": {
        "TCGA-77-7335": [
            ("chr3", 52053983, "G", "GCTC"),
            ("chr3", 52404546, "A", "ACT"),
        ],
        "TCGA-94-7557": [("chr8", 8328077, "C", "CT")],
        "TCGA-56-7582": [
            ("chr21", 42346531, "T", "TGGTACACACACCCATACATGGGTACA"),
        ],
    },
    "synthetic": {
        "synthetic": [
            ("chr8", 8328070, "G", "T"),
            ("chr8", 8328071, "C", "A"),
            ("chr8", 8328072, "C", "A"),
            ("chr3", 52404546, "AC", "GTT"),
        ],
    },
}


def copy_pass_records(src, dest):
    """Write the header and the PASS records of the VCF src to dest."""
    with open(src) as fin, open(dest, "w") as fout:
        for line in fin:
            fields = line.rstrip("\r\n").split("\t")
            if line.startswith("#") or (len(fields) > 6 and fields[6] == "PASS"):
                fout.write(line)


def write_vcf(path, records):
    with open(path, "w") as f:
        f.write(VCF_HEADER)
        for chrom, pos, ref, alt in records:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")


def run(project, genome, input_dir, out_dir):
    matGen.SigProfilerMatrixGeneratorFunc(
        project,
        genome,
        os.path.join(input_dir, ""),
        exome=False,
        plot=False,
        seqInfo=True,
        output_directory=out_dir,
    )


def main(extdata_dir, out_dir):
    inputs = {}
    for project in ["vcf", "maf", "insfix", "synthetic"]:
        inputs[project] = os.path.join(out_dir, "input_" + project)
        os.makedirs(inputs[project], exist_ok=True)
    for name in sorted(os.listdir(extdata_dir)):
        path = os.path.join(extdata_dir, name)
        if name.endswith(".vcf"):
            copy_pass_records(path, os.path.join(inputs["vcf"], name))
        elif name.endswith(".maf"):
            shutil.copy(path, os.path.join(inputs["maf"], name))
    for project, samples in EXTRA.items():
        for sample, records in samples.items():
            write_vcf(os.path.join(inputs[project], sample + ".vcf"), records)

    for project, input_dir in inputs.items():
        run(project, "GRCh38", input_dir, os.path.join(out_dir, "spmg_" + project))


if __name__ == "__main__":
    if len(sys.argv) == 6 and sys.argv[1] == "--dir":
        input_dir, project, genome, out_dir = sys.argv[2:]
        run(project, genome, input_dir, out_dir)
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        sys.exit(
            "Usage: python run_spmg.py <extdata_dir> <out_dir>\n"
            "       python run_spmg.py --dir <input_dir> <project> <genome> "
            "<out_dir>"
        )
