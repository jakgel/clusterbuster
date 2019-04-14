#!/usr/bin/env sh

conda install pandas
conda install ephem
conda install opencv
conda install -c conda-forge pypdf2
pip3 install NFW

echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
echo "export PYTHONPATH=\PYTHONPATH:$(pwd)" >> ~/.bashrc