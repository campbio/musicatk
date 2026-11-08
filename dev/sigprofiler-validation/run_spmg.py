"""Count the musicatk example variants with SigProfilerMatrixGenerator.

One-off validation for issue #218; not part of the package. See README.md
in this folder for how to set up the environment and rerun.

Usage:
    python run_spmg.py <extdata_dir> <out_dir>

<extdata_dir> is musicatk's inst/extdata. The four VCFs are copied with
only their PASS records, because SPMG ignores the FILTER column while
musicatk keeps PASS records by default; this keeps filtering out of the
comparison. The LUSC MAF is copied unchanged. SPMG runs once on the VCFs
(project "vcf") and once on the MAF (project "maf"), on GRCh38.
"""

import os
import shutil
import sys

from SigProfilerMatrixGenerator.scripts import (
    SigProfilerMatrixGeneratorFunc as matGen,
)


def copy_pass_records(src, dest):
    """Write the header and the PASS records of the VCF src to dest."""
    with open(src) as fin, open(dest, "w") as fout:
        for line in fin:
            if line.startswith("#") or line.split("\t")[6] == "PASS":
                fout.write(line)


def main(extdata_dir, out_dir):
    vcf_dir = os.path.join(out_dir, "input_vcf")
    maf_dir = os.path.join(out_dir, "input_maf")
    os.makedirs(vcf_dir, exist_ok=True)
    os.makedirs(maf_dir, exist_ok=True)
    for name in sorted(os.listdir(extdata_dir)):
        path = os.path.join(extdata_dir, name)
        if name.endswith(".vcf"):
            copy_pass_records(path, os.path.join(vcf_dir, name))
        elif name.endswith(".maf"):
            shutil.copy(path, os.path.join(maf_dir, name))

    for project, input_dir in [("vcf", vcf_dir), ("maf", maf_dir)]:
        matGen.SigProfilerMatrixGeneratorFunc(
            project,
            "GRCh38",
            input_dir + os.sep,
            exome=False,
            plot=False,
            seqInfo=True,
            output_directory=os.path.join(out_dir, "spmg_" + project),
        )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python run_spmg.py <extdata_dir> <out_dir>")
    main(sys.argv[1], sys.argv[2])
