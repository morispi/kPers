#!/bin/bash

cd minimap2
make
cd ../
mkdir RawLongReads
tar xvzf LR.fasta.tar.gz
./prepareRawLongReads.py LR.fasta RawLongReads/
./minimap2/minimap2 -Xw5 -m100 -g10000 LR.fasta LR.fasta > LR.paf
sort -k1,1 -k3,3n -k4,4n -k8,8n -k9,9n LR.paf > tmp
mv tmp LR.paf
g++ -O2 -mtune=native -march=native -std=c++14 kPers.cpp -lrt
mv a.out kPers
