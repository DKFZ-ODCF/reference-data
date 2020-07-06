__author__ = "Fritjof Lammers"
__copyright__ = "Copyright 2019, DKFZ"
__email__ = "f.lammers@dkfz-heidelberg.de"
__license__ = ""  # TODO

"""
Wrapper to add sequence (such as PhiX-spikeins) to fasta file.
snakemake.input.ref: unmodified reference fasta
snakemake.input.add: file to be added
snakemake.output: filename of modified reference fasta
"""

from snakemake.shell import shell
shell(
    "cat {snakemake.input.ref} {snakemake.input.add} > {snakemake.output}"
)