#!/usr/bin/env bash

# Download and unzip IGV browser
wget http://data.broadinstitute.org/igv/projects/downloads/IGV_2.3.74.zip
tar -xvzf IGV_2.3.74.zip
rm IGV_2.3.74.zip
mv IGV_2.3.74 igv
