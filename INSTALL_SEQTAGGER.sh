#!/bin/bash
#
# This script will download and install seqtagger
# params:
# 1 user
# 2 password

if [ x"$1" == x ]; then
        echo "please enter username"
	exit
fi
if [ x"$2" == x ]; then
        echo "please enter password"
        exit
fi

wget --user $1 --password $2 https://biocore.crg.eu/public/mop3/seqtagger/SeqTagger.tar.gz
tar -zvxf SeqTagger.tar.gz

rm SeqTagger/extract_sequence_from_fastq.py
mv SeqTagger/models mop_preprocess/seqtagger_models
mv SeqTagger/* mop_preprocess/bin/

rm SeqTagger.tar.gz
rm -fr SeqTagger
