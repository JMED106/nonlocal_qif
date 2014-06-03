#!/bin/bash

echo "Compiling..."
cd ~/Dropbox/Doctorado/Nonlocal/simulation
cd code/src
make -s
#make clean
mv nonlocal* ../../exec
cd ../../exec
cp ../code/parameters/* ./
