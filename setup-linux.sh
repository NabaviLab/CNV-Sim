#!/usr/bin/env bash

# python
sudo apt-get install python

# pip
sudo apt-get install python-pip

# wget
sudo apt-get install wget

# git
sudo apt-get install git

# Install pysam library
sudo pip install pysam

# Install numpy library
sudo pip install numpy

# Install Bio Python
pip install biopython

# install bedtools
sudo apt-get install bedtools

# install ART
mkdir lib/ART
wget -O lib/ART/art.tgz http://www.niehs.nih.gov/research/resources/assets/docs/artsrcgreatsmokymountains041716linuxtgz.tgz
tar -xvzf lib/ART/art.tgz -C lib/ART
rm lib/ART/art.tgz
mv lib/ART/art_src_GreatSmokyMountains_Linux/art_illumina lib/ART/art_illumina
rm -rf lib/ART/art_src_GreatSmokyMountains_Linux
