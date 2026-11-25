#!/bin/bash

# General

## html utilities

function make_html () {
	mkdir tmp/data/html/$1
}

## Print table header

function table_columns () {
	head -n1 $1 | sed 's/\t/\n/g' | awk '{print NR,$1}'
}

## check table header

function header_checker () {
	head $1 | column -t -s "$(printf '\011')"
}

# Retrieve genomes

## Retrieve file path of GTDB genomes from GTDB accession

function gtdbacc2path () {
	sed 's/^GB_//; s/^RS_//; s/$/_genomic.fna.gz/' | sort | join -1 1 -2 1 -a 1 - raw/gtdb/gtdb_genomes_reps_r202/genome_paths.tsv | awk '{print $2}'
}

# Retrieve proteins

## Retrieve proteins fasta files of selected genomes (GTDB accessions) and combine

function cat_proteins_bac () {
	xargs -I {} sh -c "cat gtdb_release202/protein_faa_reps/bacteria/{}_protein.faa"
}

function cat_proteins_ar () {
	xargs -I {} sh -c "cat gtdb_release202/protein_faa_reps/archaea/{}_protein.faa"
}

function retrieve_proteins_by_acc () {
	acc_tax=(`cat data/tsv/genomes_analyzed_bac.tsv data/tsv/genomes_analyzed_ar.tsv | awk -F "\t" -v OFS="\t" -v acc=$1 '$1 == acc {print $1,$111}'`);
	if [ ${acc_tax[1]} = "d__Bacteria" ]; then
		echo "${acc_tax[0]}" | xargs -I {} sh -c "cat raw/gtdb/protein_faa_reps/bacteria/{}_protein.faa"
	elif [ ${acc_tax[1]} = "d__Archaea" ]; then
		echo "${acc_tax[0]}" | xargs -I {} sh -c "cat raw/gtdb/protein_faa_reps/archaea/{}_protein.faa"
	fi
}

## Change protein names to tsv format

function protein_names2tsv () {
	sed 's/ # /\t/g; s/ID=/\t/; s/;partial=/\t/; s/;start_type=/\t/; s/;rbs_motif=/\t/; s/;rbs_spacer=/\t/; s/;gc_cont=/\t/'
}

## Make protein tsv with genomic accession, prodigal info and sequences

function proteins2tsv_bac () {
	xargs -I {} sh -c "cat raw/gtdb/protein_faa_reps/bacteria/{}_protein.faa | seqkit fx2tab | sed 's/^/{}\t/'"
}

function proteins2tsv_ar () {
	xargs -I {} sh -c "cat raw/gtdb/protein_faa_reps/archaea/{}_protein.faa | seqkit fx2tab | sed 's/^/{}\t/'"
}

function rename_nucleicacids_bac () {
	xargs -I {} sh -c "cat raw/gtdb/protein_fna_reps/bacteria/{}_protein.fna | perl -ne 's/^(>.+?) .*/\$1/;print'"
}

function rename_nucleicacids_ar () {
	xargs -I {} sh -c "cat raw/gtdb/protein_fna_reps/archaea/{}_protein.fna | perl -ne 's/^(>.+?) .*/\$1/;print'"
}

# MMSeqs

## Options for making tsv

mmseqs_tsv_opt="--format-mode 0 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,qcov,tlen,tcov"

# DIAMOND

## Options for output format

diamond_format_opt="-f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp slen scovhsp"

# List OGs assigned to query proteins

function list_cogs_hit () {
	tail -n +2 data/tsv/proteins_raw_gtdbtax_cog2020.tsv | cut -f 27 | sort | uniq
}

function list_kos_hit () {
	tail -n +2 data/tsv/proteins_raw_gtdbtax_ko.tsv | cut -f 27 | sort | uniq | sed 's/^ko://'
}

# Retrieve genomic context

function retrieve_context () {
	bash scripts/bash/retrieve_genomic_context_cog.sh $1 $2
}
