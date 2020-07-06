__author__ = "Fritjof Lammers"
__copyright__ = "Copyright 2019, DKFZ"
__email__ = "f.lammers@dkfz-heidelberg.de"
__license__ = ""  # TODO

"""
Wrapper to download using wget to a temporaryfile and gunzip the contents
Using a tmpfile is necessary because building a curl-based pipe was not possible due to difficulties
using curl on ftp links with the proxy-server
"""

from snakemake.shell import shell
shell(
    "tmpfile=$(mktemp /tmp/snakemake-wrapper-rsync-curl.XXXXXXX) && "
    "wget -e use_proxy=yes -O $tmpfile {snakemake.params.link} && " 
    "gunzip -c $tmpfile > {snakemake.output} && "
    "rm $tmpfile "
)