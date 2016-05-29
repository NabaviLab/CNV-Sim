#!/usr/bin/env bash

# Download and unzip BWA
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.13.tar.bz2
tar -xvzf bwa-0.7.13.tar.bz2
rm bwa-0.7.13.tar.bz2
mv bwa-0.7.13/ bwa

