FROM ubuntu:14.04
MAINTAINER Abdelrahman Hosny <abdelrahman.hosny@hotmail.com>

RUN apt-get update && \
    apt-get install -y build-essential python-pip wget bedtools
RUN apt-get install -y zlib1g-dev libcurl4-openssl-dev python-dev libxml2-dev libxslt-dev
RUN pip install pysam && \
    pip install numpy && \
    pip install biopython

RUN mkdir -p cnvsim/ART && \
    wget -O cnvsim/ART/art.tgz http://www.niehs.nih.gov/research/resources/assets/docs/artbingreatsmokymountains041716linux32tgz.tgz && \
    tar -xvzf cnvsim/ART/art.tgz -C cnvsim/ART && \
    rm cnvsim/ART/art.tgz && \
    mv cnvsim/ART/art_bin_GreatSmokyMountains/art_illumina cnvsim/ART/art_illumina && \
    rm -rf cnvsim/ART/art_src_GreatSmokyMountains_Linux

COPY cnvsim/Wessim /cnvsim/Wessim
COPY cnvsim/exome_simulator.py /cnvsim/exome_simulator.py
COPY cnvsim/genome_simulator.py /cnvsim/genome_simulator.py
COPY cnvsim/fileio.py /cnvsim/fileio.py
COPY cnvsim/__init__.py /cnvsim/__init__.py
COPY cnv-sim.py /