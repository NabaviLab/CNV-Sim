FROM ubuntu:14.04
MAINTAINER Abdelrahman Hosny <abdelrahman.hosny@hotmail.com>

RUN apt-get update && \
    apt-get install -y build-essential python-pip wget bedtools
RUN apt-get install -y zlib1g-dev libcurl4-openssl-dev python-dev libxml2-dev libxslt-dev
RUN pip install pysam && \
    pip install numpy && \
    pip install biopython

RUN mkdir lib/ART && \
    wget -O lib/ART/art.tgz http://www.niehs.nih.gov/research/resources/assets/docs/artbingreatsmokymountains041716linux32tgz.tgz && \
    tar -xvzf lib/ART/art.tgz -C lib/ART && \
    rm lib/ART/art.tgz && \
    mv lib/ART/art_bin_GreatSmokyMountains/art_illumina lib/ART/art_illumina && \
    rm -rf lib/ART/art_src_GreatSmokyMountains_Linux

COPY lib/Wessim /lib/Wessim
COPY lib/exome_simulator.py /lib/exome_simulator.py
COPY lib/genome_simulator.py /lib/genome_simulator.py
COPY lib/fileio.py /lib/fileio.py
COPY lib/__init__.py /lib/__init__.py
COPY cnv-sim.py /