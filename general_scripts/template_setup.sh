#!/bin/bash

echo "Input directory where template needs to be created..."

read path

echo "Input content for the readme file..."

read readme

echo "Input content for the template file..."

read commands

cd

cd "$path"

mkdir project
mkdir project/scripts
mkdir project/bash
echo '#!/bin/bash'> project/bash/runall.sh
mkdir project/input
mkdir project/output
mkdir project/logs
echo $readme> project/readme.txt

mkdir datatsets

mkdir general_scripts

mkdir bin


