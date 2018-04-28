#!/bin/bash

echo $1 $2

python ConvertBlasrToPaf.py "$1" "$2"
./sortPaf.sh "$2" "TMPSORTEDPAF"
mv "TMPSORTEDPAF" "$2"
