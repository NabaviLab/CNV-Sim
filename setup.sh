#!/usr/bin/env bash

mkdir tmp
cd tmp

# INSTALL GENERAL DEPENDENCIES (Mac OS)

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

# java
brew tap caskroom/cask
brew install brew-cask
brew cask install java