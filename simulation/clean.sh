#!/bin/bash
cd ~/Dropbox/Doctorado/Nonlocal/simulation

echo "Removing all data and executables from $(pwd)/exec"
read -p "Are you sure?" -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm -r $(pwd)/exec/*
else
    exit 1
fi

echo "Done."
exit 0