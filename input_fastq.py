#!/bin/env python

"""
Functions to parse directories of input fastq files, create symlinks with expected filename structures, and return metadata.

Clare Sloggett, VLSCI
"""

import sys
import re
import os.path
from collections import defaultdict
from rubra.utils import (mkLink)

def parse_and_link(file, symlink_dir, metadata_dict):
    """
    Parse metadata out of input filename and construct symlink.
    Takes a fastq filename, destination directory, and a metadata dict, which should be of type defaultdict(dict).
    Parse the filename to get information on the sample name, run, read #, etc.
    Medadata is added to the provided metadata_dict.
    Some metadata is used to build symlinks, to guarantee filename uniqueness and a regular naming structure.\
    Currently parsing by assuming AGRF naming structure and paired-end reads
    Currently will ONLY handle gzipped files, to avoid multiple links to the same data.
    """
    match_old = re.match(r".*?/([^_/]+)_([a-zA-Z0-9-.]+)_s_([0-9]+)_(1|2)_sequence.txt.gz",file)
    match_new = re.match(r".*?/([a-zA-Z0-9-.]+)_([^_/]+)_[CAGTN]+_L([0-9]+)_R(1|2).fastq.gz",file)
    if match_old:
        run_id = match_old.group(1)
        sample = match_old.group(2)
        lane = int(match_old.group(3))
        pair = match_old.group(4)
        encoding = 'I'
    elif match_new:
        run_id = match_new.group(2)
        sample = match_new.group(1)
        lane = int(match_new.group(3))
        pair = match_new.group(4)
        encoding = 'S'
    else:
        print "Unable to parse name of fastq file %s ." % file
        sys.exit(1)
    newfile = os.path.join(symlink_dir, "%s_%s_L%d_%s.fastq.gz" % 
                                    (sample, run_id, lane, pair))
    metadata_dict[os.path.basename(newfile)]['sample'] = sample
    metadata_dict[os.path.basename(newfile)]['run_id'] = run_id
    metadata_dict[os.path.basename(newfile)]['lane'] = lane
    metadata_dict[os.path.basename(newfile)]['pair'] = pair
    metadata_dict[os.path.basename(newfile)]['encoding'] = encoding
    relative_sourcefile = os.path.relpath(file, symlink_dir)
    mkLink(relative_sourcefile, newfile)
    return newfile

