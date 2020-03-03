#! /usr/bin/bash

# get working directory
WORK_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

genome_DIR=/mnt/imaging.data/jsemple/genomeVer/ws265
export genomeIdx=${genome_DIR}/ce11_genome.fa
blacklist=${genome_DIR}/ce11-blacklist.v2.bed


