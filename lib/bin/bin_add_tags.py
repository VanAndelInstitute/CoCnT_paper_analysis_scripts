#!/usr/bin/env python3

import argparse
from sys import stdout
from random import choice
from simplesam import Reader, Writer

def parse_qname(ss, index=None):
    substring = ss.qname.split(':')[4]
    substring = '_'.join(substring.split('_')[1:5])
    if index is not None:
        substring = f"{substring}-{index}"
    return {'CB': substring, 'CR': substring }

def iterate(args):
    with Reader(args.bam) as bam, Writer(stdout, bam.header) as stdout_sam:
        bam.header.get('@HD')['VN:1.4'] = ['SO:coordinate']
        for read in bam:
            read.tags.update(parse_qname(read, index=args.index))
            stdout_sam.write(read)

def main():
    parser = argparse.ArgumentParser(prog='addTags', description="parse BAM sequence name for barcode and add as bam tags")
    ## arguments

    ## bam file
    parser.add_argument('bam', type=argparse.FileType('r'), help=" BAM file ")

    ## index appended to barcode
    parser.add_argument('-i', '--index', required=False, help="index appended to barcode")

    parser.set_defaults(func=iterate)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
