import sys

fin = open(sys.argv[1])
fout = open(sys.argv[2], 'w')

l = fin.readline()
while l != '':
	t = l.split(" ")
	qName = t[0].split("/")[0]
	tName = t[1]
	qStrand = t[2]
	tStrand = t[3]
	score = t[4]
	percentSimilarity = t[5]
	tStart = t[6]
	tEnd = t[7]
	tLength = t[8]
	qStart = t[9]
	qEnd = t[10]
	qLength = t[11]
	nCells = t[12]
	if qStrand == tStrand:
		strand = '+'
	else:
		strand = '-'
	if qName != tName:
		fout.write(qName + "\t" + qLength + "\t" + qStart + "\t" + qEnd + "\t" + strand + "\t" + tName + "\t" + tLength + "\t" + tStart + "\t" + tEnd + "\t" + "0" + "\t" + "0" + "\t" + "0" + "\n")
	l = fin.readline()
