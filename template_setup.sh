#!/bin/bash

echo "Input directory where template needs to be created..."

read path

echo "Input content for the readme file..."

read readme

echo "Input content for the template file..."

read commands

cd

cd $path

mkdir projects
mkdir projects/scripts
mkdir projects/bash
echo '#!/bin/bash'> projects/bash/runall.sh
mkdir projects/input
mkdir projects/output
echo $readme> projects/readme.txt
echo $commands> projects/commands.txt

mkdir datatsets

mkdir general_scripts

mkdir bin


