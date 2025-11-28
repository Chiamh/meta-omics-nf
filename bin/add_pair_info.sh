#!/bin/bash
set -e
set -u

read_1=${1}
read_2=${2}

sed "s/\s\([0-9]\+\):N:/\/1/g" -i ${read_1} && gzip ${read_1}
sed "s/\s\([0-9]\+\):N:/\/2/g" -i ${read_2} && gzip ${read_2}
