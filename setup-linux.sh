#!/usr/bin/env bash

# python
sudo apt-get install python

# pip
sudo apt-get install python-pip

# wget
sudo apt-get install wget

# Install pysam library
sudo pip install pysam

# Install numpy library
sudo pip install numpy

# Install Bio Python
sudo pip install biopython

# install bedtools
sudo apt-get install bedtools

# install ART
mkdir lib/ART
wget -O cnvsim/ART/art.tgz http://www.niehs.nih.gov/research/resources/assets/docs/artsrcgreatsmokymountains041716linuxtgz.tgz
tar -xvzf cnvsim/ART/art.tgz -C cnvsim/ART
rm cnvsim/ART/art.tgz
mv cnvsim/ART/art_src_GreatSmokyMountains_Linux/art_illumina cnvsim/ART/art_illumina
rm -rf cnvsim/ART/art_src_GreatSmokyMountains_Linux
