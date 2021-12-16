###	ybyra v0.2
###
### simple exploratory Y-chromosome caller
###
### TP 06/21 (thomaz.pinotti AT gmail.com)

BAM = [line.rstrip() for line in open("ysamples.bam")] # list bam files
SITES = "/path/to/SNPindex/10Mb_callable/SNPindex_all-10Mb.angsd" # change as desired
MIN_DEPTH = "1" # minimum read depth at a region to call a variant
TRIM = "0" # number of bases to trim at read termini
MIN_MAPQ = "30" # minimum mapQ
MIN_BQ = "20" # minimum baseQ
ISOGG_NAME = "/path/to/SNPindex/SNP-names-hap.pos"

rule all:
	input:
		"no_hits.list"

rule angsd:
	input:
		"/path/to/{bam}.bam"
	output:
		mafs = temp("{bam}.mafs.gz"),
		arg = temp("{bam}.arg"),
	params:
		sites = SITES,
		mindepth = MIN_DEPTH,
		trim = TRIM,
		minmapq = MIN_MAPQ,
		minbq = MIN_BQ,
		outname = "{bam}",
	shell:
		"angsd -sites {params.sites} -DoMajorMinor 3 -doMaf 8 -doCounts 1 -trim {params.trim} -setMinDepthInd {params.mindepth} -minMapQ {params.minmapq} -minQ {params.minbq} -out {params.outname} -i {input} -r Y: "

rule filter:
	input:
		"{bam}.mafs.gz"
	output:
		temp("{bam}.pos")
	shell:
		"""zless -S {input} | awk '$5>0.9' | cut -f 2 > {output}"""

rule panel_match:
	input:
		"{bam}.pos"
	output:
		"{bam}.hits"
	params:
		isogg = ISOGG_NAME,
	shell:
		"""
		awk 'NR==FNR{{a[$1];next}} ($1) in a' {input} {params.isogg} | 
		awk '{{print ($(NF-2) " "$0)}}'| 
		sort -k 1,1 -k 2,2 -u > {output}
		"""

rule no_hits:
	input:
		expand("{bam}.hits", bam=BAM)
	output:
		"no_hits.list"
	shell:
		"""echo "List of bam files with 0 hits in used SNP index/filters:" > no_hits.list ;"""
		"find *hits -size 0 -print >> no_hits.list ;"
		"find *hits -size 0 -delete"