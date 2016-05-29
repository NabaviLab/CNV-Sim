#!/usr/bin/env bash

mkdir tmp
cd tmp

# Install general dependencies (Mac OS only)
# Homebrew
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# Python
brew install python

# PIP
curl -O http://python-distribute.org/distribute_setup.py
python distribute_setup.py
curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
python get-pip.py

# wget
homebrew install wget

# git
brew install git

# Clean after installation
rm -t tmp
