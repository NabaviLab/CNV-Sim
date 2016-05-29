#!/usr/bin/env bash

# Download and unzip samtools
wget https://sourceforge.net/projects/samtools/files/samtools/1.3.1/samtools-1.3.1.tar.bz2
tar -xvzf samtools-1.3.1.tar.bz2
rm samtools-1.3.1.tar.bz2
mv samtools-1.3.1/ samtools

# Install
cd samtools
./configure
make
make install