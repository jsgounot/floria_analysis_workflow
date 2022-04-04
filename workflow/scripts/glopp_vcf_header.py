# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2022-04-02 11:26:26
# @Last Modified by:   jsgounot
# @Last Modified time: 2022-04-02 19:35:42

# A copy of glopp script for a better embedding 
# https://github.com/bluenote-1577/glopp/blob/flow/scripts/write_contig_headers_vcf.py


from sys import argv

def main(vcf_file, new_vcf):
	refs = set()

	with open(vcf_file) as f :

		refs = {
			line.split()[0] for line in f
			if line != "" and line[0] != '#'
		}

		refs = sorted(refs)

	with open(vcf_file) as f :
		with open(new_vcf, 'w') as fo:
			for idx, line in enumerate(f):
				if idx == 2:
					for ref in refs:
						fo.write(f"##contig=<ID={ref}>\n")

				fo.write(line)

vcf_file = snakemake.input[0]
new_vcf = snakemake.output[0]
main(vcf_file, new_vcf)