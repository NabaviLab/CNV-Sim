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

# git
brew install git

# clean after installation
rm -t tmp

# Install pysam library
pip install pysam

# Install numpy library
pip install numpy

# install bedtools
brew install homebrew/science/bedtools
