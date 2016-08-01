#!/usr/bin/env bash

# Homebrew
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# python
brew install python

# pip
curl -O http://python-distribute.org/distribute_setup.py
python distribute_setup.py
curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
python get-pip.py

# wget
brew install wget

# clean after installation
rm -t tmp

# Install pysam library
pip install pysam

# Install numpy library
pip install numpy

# Install Bio Python
pip install biopython

# install bedtools
brew install homebrew/science/bedtools

# install ART
mkdir cnvsim/ART
wget -O cnvsim/ART/art.tgz http://www.niehs.nih.gov/research/resources/assets/docs/artsrcgreatsmokymountains041716macostgz.tgz
tar -xvzf cnvsim/ART/art.tgz -C lib/ART
rm cnvsim/ART/art.tgz
mv cnvsim/ART/art_src_GreatSmokyMountains_MacOS/art_illumina cnvsim/ART/art_illumina
rm -rf cnvsim/ART/art_src_GreatSmokyMountains_MacOS


