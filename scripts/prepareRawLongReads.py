#!/usr/bin/python3

import os
import sys
import subprocess

#os.makedirs(os.path.dirname(sys.argv[1])) #, exist_ok=True)
f = open(sys.argv[1])

finalString = ""
id = f.readline()
while id != "":
	seq = f.readline()
	os.makedirs(os.path.dirname(sys.argv[2] + id[1:-1]), exist_ok=True)
	out = open(sys.argv[2] + id[1:-1], 'w')
	out.write(id + seq.lower())
	out.close()
	id = f.readline()
f.close()
