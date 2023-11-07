#!/bin/bash

echo "Check dependency"

if type zenodo_get &>/dev/null; then 
    echo "zenodo_get found"
    mkdir data
    zenodo_get 10.5281/zenodo.10036117 -o ./data

else
    echo "zenodo_get not found"
fi