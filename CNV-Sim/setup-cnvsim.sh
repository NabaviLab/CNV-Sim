#!/usr/bin/env bash

# Install pysam library
pip install pysam

# Install numpy library
pip install numpy

# Clone and organize Wessim
git clone https://github.com/sak042/Wessim.git
mv Wessim/ tmp
mv tmp/Wessim_ver_1.0/ .
mv Wessim_ver_1.0/ Wessim
rm -r -f tmp/

# install bedtools
brew install homebrew/science/bedtools
