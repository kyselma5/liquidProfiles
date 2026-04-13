#!/bin/bash

set -e

cd src

echo "BUILDING"
make clean
make

echo "RUNNING"
time ./build/out config.cfg

echo "GENERATING GRAPHS (this takes a few seconds)"
jupyter nbconvert --to notebook --execute other/graphsGenerator.ipynb --output executed.ipynb
#jupyter notebook other/executed.ipynb
echo "DONE"
