
"""
Read all output flagstat files and store values in a text table.
Creates two tables in the output directory: readcounts.txt and readcounts_fractions.txt, which contains useful human-readable percentage statistics.
Each row will be a sample.
Values to store (columns) are:
    Total reads from SAMPLE.bam.flagstat
    Mapped reads from SAMPLE.bam.flagstat
    Mapped reads from SAMPLE.dedup.bam.flagstat
There is an assumption that the bams were mapped with something like bwa, which includes unmapped reads.

Usage: python count_flagstat_wgs.py flagstat_directory output_directory
"""

import sys
import os
import re
import optparse
from collections import defaultdict

class FlagstatParseException (Exception):
    pass

def read_flagstat(filename):
    """
    Given a filename, parse the flagstat values and return as a hash.
    Relying on flagstat contents being in usual order.
    """
    numbers = re.compile(r'^(\d+)\s+\+\s+(\d+)\s+')
    values = {}
    f = open(filename)
    for field in ['total',
                    'duplicates',
                    'mapped',
                    'paired',
                    'read1',
                    'read2',
                    'properly_paired',
                    'both_mapped',
                    'singletons',
                    'mate_distant',
                    'mate_distant_goodqual']:
        line = f.readline().strip()
        match = numbers.match(line)
        if not match:
            raise FlagstatParseException
        values[field] = int(match.group(1))
        values[field+'_QCfailed'] = int(match.group(2))
    f.close()
    return values

# -----

# Get arguments and input filenames
parser = optparse.OptionParser(usage=__doc__)
#parser.add_option()
(options, args) = parser.parse_args()
if len(args) != 2:
    parser.error("Wrong number of arguments - see usage info")
in_dir = args[0]
out_dir = args[1]
if not (os.path.exists(out_dir) and os.path.isdir(out_dir)):
    sys.exit("There does not seem to be a directory %s , exiting" % out_dir)

filenames = os.listdir(in_dir)

#print ', '.join(filenames)

alignedname = re.compile('^([^_\.]+).bam.flagstat')
dedupname = re.compile('^([^_\.]+).dedup.bam.flagstat')

samples = defaultdict(dict)
for filename in filenames:
    if dedupname.match(filename):
        match = dedupname.match(filename)
        name = match.group(1)
        values = read_flagstat( os.path.join(in_dir, filename) )
        samples[name]['deduped'] = values['mapped']    
    elif alignedname.match(filename):
        match = alignedname.match(filename)
        name = match.group(1)
        values = read_flagstat( os.path.join(in_dir, filename) )
        samples[name]['mapped'] = values['mapped']
        samples[name]['total'] = values['total']
        
#print ', '.join(samples.keys())

tablefile = os.path.join(out_dir, "readcounts.txt")
tablefile_plus = os.path.join(out_dir, "readcounts_fractions.txt")

OUT_TABLE = open(tablefile, 'w')
OUT_TABLEPLUS = open(tablefile_plus, 'w')

OUT_TABLE.write("Sample\tTotal\tMapped\tDeduped\n") 
OUT_TABLEPLUS.write("Sample\tTotal\tMapped\t%\tDeduped\t%\n")

for sample in sorted(samples.keys()):
    values = samples[sample]
    fraction_mapped = float(values['mapped'])/values['total']
    fraction_deduped = float(values['deduped'])/values['mapped']

    OUT_TABLE.write( "%s\t%d\t%d\t%d\n" % (sample, values['total'], values['mapped'], values['deduped']) )
    
    OUT_TABLEPLUS.write( "%s\t%d\t%d\t%.3f\t%d\t%.3f\n" % (sample, values['total'], values['mapped'], fraction_mapped, values['deduped'], fraction_deduped) )

OUT_TABLE.close()
OUT_TABLEPLUS.close()
