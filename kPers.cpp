#include <iostream>
#include <chrono>
#include <mutex>
#include <future>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm> 
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <algorithm>
#include <string>
#include <seqan/index.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/arg_parse.h>

using namespace seqan;

// DnaString genome;
// typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;
// Index<DnaString, BidirectionalIndex<FMIndex<>>> ind;
// int maxErr;
// int nbRead;
// int lenRead;

// void loadIndex(std::set<std::string> kMers, int merSize) {
// 	// Defining constants from available data
// 	maxErr  = 2;
// 	nbRead  = kMers.size();
// 	lenRead = merSize; // supposing all sequences have the same length

// 	// Creating the sequence to index
// 	genome = "";
// 	for (std::string k : kMers) {
// 		append(genome, k);
// 	}

// 	// Setting up the index.
// 	ind = Index<DnaString, BidirectionalIndex<FMIndex<>>>(genome);
// }

struct alignment_t {
	std::string qName;
	int qLength;
	int qStart;
	int qEnd;
	bool strand;
	std::string tName;
	int tLength;
	int tStart;
	int tEnd;
	std::string resMatches;
	std::string alBlockLen;
	std::string mapQual;

	bool operator<(const alignment_t& a2) const {
		if (qName < a2.qName) {
			return true;
 		} else if (qName == a2.qName && qLength < a2.qLength) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart < a2.qStart) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd < a2.qEnd) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand < a2.strand) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName < a2.tName) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength < a2.tLength) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart < a2.tStart) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart == a2.tStart && tEnd < a2.tEnd) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart == a2.tStart && tEnd == a2.tEnd && resMatches < a2.resMatches) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart == a2.tStart && tEnd == a2.tEnd && resMatches == a2.resMatches && alBlockLen < a2.alBlockLen) {
 			return true;
 		} else if (qName == a2.qName && qLength == a2.qLength && qStart == a2.qStart && qEnd == a2.qEnd && strand == a2.strand && tName == a2.tName && tLength == a2.tLength && tStart == a2.tStart && tEnd == a2.tEnd && resMatches == a2.resMatches && alBlockLen == a2.alBlockLen && mapQual < a2.mapQual) {
 			return true;
 		} else {
 			return false;
 		}
	}
};

struct consensus_t {
	std::string readId;
	int startPos;
	bool strand;
	int relPos;

	bool operator<(const consensus_t& c2) const {
		if (readId < c2.readId) {
			return true;
		} else if (readId == c2.readId && startPos < c2.startPos) {
			return true;
		} else if (readId == c2.readId && startPos == c2.startPos && strand < c2.strand) {
			return true;
		} else if (readId == c2.readId && startPos < c2.startPos && strand == c2.strand && relPos < c2.relPos) {
			return true;
		} else {
			return false;
		}
	}
};

std::string reverseComplement(std::string seq) {
	std::string res = std::string(seq);
	for (int i = 0 ; i < seq.length() ; i++) {
		switch(seq[i]) {
			case 'A':
				res[seq.length() - i - 1] = 'T';
				break;
			case 'C':
				res[seq.length() - i - 1] = 'G';
				break;
			case 'G':
				res[seq.length() - i - 1] = 'C';
				break;
			case 'T':
				res[seq.length() - i - 1] = 'A';
				break;
			case 'a':
				res[seq.length() - i - 1] = 't';
				break;
			case 'c':
				res[seq.length() - i - 1] = 'g';
				break;
			case 'g':
				res[seq.length() - i - 1] = 'c';
				break;
			case 't':
				res[seq.length() - i - 1] = 'a';
				break;
		}
	}

	return res;
}	

std::string getCannonicalSeq(std::string seq) {
	std::string revComp = reverseComplement(seq);
	return seq < revComp ? seq : revComp;
}

std::map<std::string, int> getKPersCounts(std::set<std::string> kMers, std::set<std::string> candidates, int maxError, int merSize, std::map<std::string, int>& merCounts) {
	std::map<std::string, int> kPers;

	const int maxErr  = 2;
    const int nbRead  = candidates.size();
    // assuming all sequences have the same length
    const int lenRead = merSize;

    // Creating the sequence to index
    DnaString genome;
    for (std::string c : candidates){
        append(genome,c);
    }
    

    // Setting up the index.
    typedef FastFMIndexConfig<void, uint32_t, 2, 1> TFastConfig;
    Index<DnaString, BidirectionalIndex<FMIndex<>>> index(genome);

    int nbkpers = 0;
    int count;
    int cptpers;
    for (std::string k : kMers) {
    	count = 0;
    	cptpers = 0;
	    // Kmer to look for
	    DnaString kmer = k;

	    // hit storrage
		typedef std::map< Infix<String<Dna> >::Type , std::vector<int> >   TKmerCounter;
		TKmerCounter hits;

		// Used to store small sequences by reference to the "genome"
		Infix<DnaString>::Type inf;


		// Define what to do when a hit is found.
		auto delegate = [& hits, & lenRead, & nbRead, & genome, & inf](auto & iter, DnaString const & needle, uint8_t errors) {
		    // Storing each kmer and position  in a hit table
		    for (auto occ : getOccurrences(iter)){
		        inf = infix(genome,occ,occ+16);
		        std::vector<int> currentVector = hits[inf];
		        if( occ%lenRead == 0 and  std::find(currentVector.begin(), currentVector.end(), occ) == currentVector.end() ){
		            currentVector.push_back(occ);
		        }


		        hits[inf] = currentVector;
		    }
		};	

	    // Actually performing the research
	    find<0, maxErr>(delegate, index, kmer, EditDistance());
	    
	    // Processing results
	    for(auto elem: hits){
	        if(length(elem.second) != 0)   {
	        	String<char> char_target;
	        	assign(char_target, elem.first);
	        	std::string s;
	        	for (char c : char_target) {
	        		s += c;
	        	}
	        	std::transform(s.begin(), s.end(), s.begin(), ::tolower);
	        	count += merCounts[s];
	        	cptpers++;
	        }
	    }
	    kPers[k] = count;
	    nbkpers++;
	    // std::cerr << k << " has " << cptpers << " k-pers" << " and count " << count << std::endl;
	    // std::cerr << "computed " << nbkpers << " k-pers out of " << kMers.size() << std::endl;
	}

	return kPers;
}

std::map<std::string, int> getKMersCounts(std::set<std::string> sequences, int merSize) {
	std::map<std::string, int> merCounts;
	int i;
	std::string m;

	for (std::string s : sequences) {
		i = 0;
		while (s.length() >= merSize && i < s.length() - merSize + 1) {
			m = getCannonicalSeq(s.substr(i, merSize));
			merCounts[m]++;
			i++;
		}
	}

	return merCounts;
}

std::pair<std::set<std::string>, std::set<std::string>> computeGenomicKMers(std::set<std::string> similarRegions, int merSize, int thFreq, int maxError, int thPers) {
	std::set<std::string> genMers;
	std::set<std::string> errMers;
	std::set<std::string> candidates;
	std::map<std::string, int> persCounts;

	// k-mers computation
	std::map<std::string, int> merCounts = getKMersCounts(similarRegions, merSize);	
	for (std::map<std::string, int>::iterator it = merCounts.begin(); it != merCounts.end(); it++) {
		if (it->second >= thFreq) {
			genMers.insert(it->first);
		} else if (it->second > 1) {
			errMers.insert(it->first);
		}
		candidates.insert(it->first);
	}

	// std::cerr << genMers.size() << " genomic k-mers" << std::endl;
	// std::cerr << errMers.size() << " k-pers to compute" << std::endl;
	// std::cerr << candidates.size() << " candidates for k-pers" << std::endl;

	// k-pers computation
	persCounts = getKPersCounts(errMers, candidates, maxError, merSize, merCounts);
	int addCpt = 0;
	for (std::map<std::string, int>::iterator it = persCounts.begin(); it != persCounts.end(); it++) {
		if (it->second >= thPers) {
			addCpt++;
			genMers.insert(it->first);
		} 
		// else {
		// 	errMers.insert(it->first);
		// }
	}
	return std::make_pair(genMers, errMers);
}

std::set<std::string> processRegion(std::set<consensus_t> alignments, int SHIFT, int length, std::string readsDir) {
	std::set<std::string> sequences;
	std::string curHeader, line, region;
	int beg, end;
	
	for (std::set<consensus_t>::iterator it = alignments.begin(); it != alignments.end(); it++) {
		std::ifstream f(readsDir + it->readId);
		if (!f.fail()) {
			getline(f, line);
			getline(f, line);
			f.close();
			if (it->strand == true) {
				line = reverseComplement(line);
			}
			beg = it->startPos + SHIFT - it->relPos;
			if (beg < line.length()) {
				region = line.substr(beg, length);
				sequences.insert(region);
			}
		}
	}

	return sequences;
}

consensus_t createConsensusStruct(std::string name, int qStart, bool strand, int posrel) {
	consensus_t t;
	t.readId = name;
	t.startPos = qStart;
	t.strand = strand;
	t.relPos = posrel;

	return t;
}

std::set<std::string> processRead(std::set<alignment_t> alignments, std::string readsDir, int merSize, int minSupport, int minRegionLength, int thFreq, int maxError, int thPers) {
	std::set<std::string> genMers;
	std::set<std::string> errMers;
	std::set<std::string> regSeqs;
	std::pair<std::set<std::string>, std::set<std::string>> allMers;
	std::set<consensus_t> curRegion;
	int min, max, sup, startPoint;
	min = -1;
	max = -1;
	sup = 0;

	for (std::set<alignment_t>::iterator it = alignments.begin(); it != alignments.end(); it++) {
		// std::cerr << "processing an alignment" << std::endl;
		if (min == -1 and max == -1) {
			// std::cerr << "1" << std::endl;
			startPoint = it->qStart;
			min = it->qStart;
			max = it->qEnd;
			curRegion.insert(createConsensusStruct(it->qName, it->qStart, false, it->qStart));
			curRegion.insert(createConsensusStruct(it->tName, it->tStart, it->strand, it->qStart));
			sup += 2;
			//elif (int(al[2]) > min and int(al[2]) > max - minRegLen) or        (int(al[3]) < max and min > int(al[3]) - minRegLen) or (int(al[2]) > min and int(al[3]) < max and int(al[2]) > int(al[3]) - minRegLen):
		} else if ((it->qStart > min and it->qStart > max - minRegionLength) or (it->qEnd < max and min > it->qEnd - minRegionLength) or (it->qStart > min and it->qEnd < max and it->qStart > it->qEnd - minRegionLength)) {
			// std::cerr << "2" << std::endl;
			if (sup >= minSupport) {
				regSeqs = processRegion(curRegion, min, max - min, readsDir);
				allMers = computeGenomicKMers(regSeqs, merSize, thFreq, maxError, thPers);
				genMers.insert(allMers.first.begin(), allMers.first.end());
				errMers.insert(allMers.second.begin(), allMers.second.end());
			}
			startPoint = it->qStart;
			min = it->qStart;
			max = it->qEnd;
			sup = 2;
			curRegion.clear();
			curRegion.insert(createConsensusStruct(it->qName, it->qStart, false, it->qStart));
			curRegion.insert(createConsensusStruct(it->tName, it->tStart, it->strand, it->qStart));
		} else if (it->qStart <= min) {
			// std::cerr << "3" << std::endl;
			if (it->qEnd < max) {
				if (sup >= minSupport and max - it->qEnd >= minRegionLength) {
					regSeqs = processRegion(curRegion, it->qEnd, max - it->qEnd, readsDir);
					allMers = computeGenomicKMers(regSeqs, merSize, thFreq, maxError, thPers);
					genMers.insert(allMers.first.begin(), allMers.first.end());
					errMers.insert(allMers.second.begin(), allMers.second.end());
				}
				max = it->qEnd;
			}
			curRegion.insert(createConsensusStruct(it->tName, it->tStart, it->strand, it->qStart));
			sup += 1;
		} else if (it->qStart < max) {
			// std::cerr << "4" << std::endl;
			if (sup >= minSupport and it->qStart - min >= minRegionLength) {
				regSeqs = processRegion(curRegion, min, it->qStart - min, readsDir);
				allMers = computeGenomicKMers(regSeqs, merSize, thFreq, maxError, thPers);
				genMers.insert(allMers.first.begin(), allMers.first.end());
				errMers.insert(allMers.second.begin(), allMers.second.end());
			}
			min = it->qStart;
			if (it->qEnd < max) {
				if (sup >= minSupport and max - it->qEnd >= minRegionLength) {
					regSeqs = processRegion(curRegion, it->qEnd, max - it->qEnd, readsDir);
					allMers = computeGenomicKMers(regSeqs, merSize, thFreq, maxError, thPers);
					genMers.insert(allMers.first.begin(), allMers.first.end());
					errMers.insert(allMers.second.begin(), allMers.second.end());
				}
				max = it->qEnd;
			}
			curRegion.insert(createConsensusStruct(it->tName, it->tStart, it->strand, it->qStart));
			sup += 1;
		}
		// curRegion.insert(createConsensusStruct(it->tName, it->tStart, it->strand, it->qStart));
		// sup += 1;
	}
	// std::cerr << "out of the loop" << std::endl;
	if (sup >= minSupport) {
		regSeqs = processRegion(curRegion, min, max - min, readsDir);
		allMers = computeGenomicKMers(regSeqs, merSize, thFreq, maxError, thPers);
		genMers.insert(allMers.first.begin(), allMers.first.end());
		errMers.insert(allMers.second.begin(), allMers.second.end());
	}

	return genMers;
}

alignment_t getAlignmentFromString(std::string al) {
	alignment_t t;
	std::string token;
	std::stringstream iss(al);
	getline(iss, token, '\t');
	t.qName = token;
	getline(iss, token, '\t');
	t.qLength = stoi(token);
	getline(iss, token, '\t');
	t.qStart = stoi(token);
	getline(iss, token, '\t');
	t.qEnd = stoi(token);
	getline(iss, token, '\t');
	t.strand = token == "+" ? false : true;
	getline(iss, token, '\t');
	t.tName = token;
	getline(iss, token, '\t');
	t.tLength = stoi(token);
	getline(iss, token, '\t');
	t.tStart = stoi(token);
	getline(iss, token, '\t');
	t.tEnd = stoi(token);
	getline(iss, token, '\t');
	t.resMatches = token;
	getline(iss, token, '\t');
	t.alBlockLen = token;
	getline(iss, token, '\t');
	t.mapQual = token;

	return t;
}

void goodRegionsGenK(std::string alignmentFile, std::string readsDir, int merSize, int minSupport, int minRegionLength, int thFreq, int maxError, int nbThreads, int thPers) {
	// loadIndex();
	std::ifstream f(alignmentFile);
	std::set<alignment_t> curReadAlignments;
	alignment_t al;
	std::set<std::string> genMers;
	std::string curRead, line;
	curRead = "";
	int id = 0;
	int readsNumber = 0;
	int i = 0;

	getline(f, line);
	while(line.length() > 0 or !curReadAlignments.empty()) {
		if (line.length() > 0) {
			al = getAlignmentFromString(line);
		}
		if (line.length() > 0 and (curRead == "" or al.qName == curRead)) {
			i++;
			curRead = al.qName;
			curReadAlignments.insert(al);
			getline(f, line);
		} else {
			genMers = processRead(curReadAlignments, readsDir, merSize, minSupport, minRegionLength, thFreq, maxError, thPers);
			std::cerr << genMers.size() << std::endl;
			for (std::string s : genMers) {
				std::cout << ">kmer" << std::endl << s << std::endl;
			}			
			curReadAlignments.clear();
			curRead = "";
		}
	}
}

std::set<std::string> processCluster(std::set<std::string> clusterIds, std::string readsDir, int merSize, int thFreq, int maxError, int thPers) {
	std::pair<std::set<std::string>, std::set<std::string>> allMers;
	std::set<std::string> clusterSeqs;
	std::string line;
	for (std::string id : clusterIds) {
		std::ifstream f(readsDir + id);
		std::pair<std::set<std::string>, std::set<std::string>> allMers;
		if (!f.fail()) {
			getline(f, line);
			getline(f, line);
			clusterSeqs.insert(line);
		}
	}

	allMers = computeGenomicKMers(clusterSeqs, merSize, thFreq, maxError, thPers);
	return allMers.first;
}

void clusterGenK(std::string alignmentFile, std::string readsDir, int merSize, int thFreq, int maxError, int nbTheads, int thPers) {
	std::ifstream f(alignmentFile);
	std::set<std::string> cluster;
	alignment_t al;
	std::set<std::string> genMers;
	std::string curRead, line;
	curRead = "";
	int id = 0;
	int readsNumber = 0;
	int i = 0;
	getline(f, line);

	while (line.length() > 0 or !cluster.empty()) {
		if (line.length() > 0) {
			al = getAlignmentFromString(line);
		}
		if (line.length() > 0 and (curRead == "" or al.qName == curRead)) {
			i++;
			curRead = al.qName;
			if (cluster.empty()) {
				cluster.insert(al.qName);
			}
			cluster.insert(al.tName);
			getline(f, line);
		} else {
			genMers = processCluster(cluster, readsDir, merSize, thFreq, maxError, thPers);
			std::cerr << genMers.size() << std::endl;
			for (std::string s : genMers) {
				std::cout << ">kmer" << std::endl << s << std::endl;
			}
			cluster.clear();
			curRead = "";
		}
	}
}

int main(int argc, char* argv[]) {
	if (argc < 2) {
		fprintf(stderr, "Usage: %s [-a alignmentFile.paf] [-d RawLongReadsDir] [-k merSize] [-s minSupportForGoodRegions] [-l minLengthForGoodRegions] [-f freqThresholdForKMers] [-e maxError] [-p freqThresholdForKPers] [-m mode (0 for regions, 1 for cluster)] [-j threadsNb] \n\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	std::string alignmentFile, readsDir;
	int merSize, minSupport, minRegionLength, thFreq, maxError, nbThreads, thPers, opt, mode;

	readsDir =  "RawLongReads/";
	minSupport = 3;
	minRegionLength = 100;
	thFreq = 3;
	maxError = 2;
	nbThreads = 1;
	thPers = 5;


	while ((opt = getopt(argc, argv, "a:d:k:s:l:f:e:p:m:j")) != -1) {
        switch (opt) {
			case 'a':
				alignmentFile = optarg;
				break;
			case 'd':
				readsDir = optarg;
				break;
			case 'k':
				merSize = atoi(optarg);
				break;
			case 's':
				minSupport = atoi(optarg);
				break;
			case 'l':
				minRegionLength = atoi(optarg);
				break;
			case 'f':
				thFreq = atoi(optarg);
				break;
			case 'e':
				maxError = atoi(optarg);
				break;
			case 'p':
				thPers = atoi(optarg);
				break;
			case 'j':
				nbThreads = atoi(optarg);
				break;
			case 'm':
				mode = atoi(optarg);
				break;
			default: /* '?' */
				fprintf(stderr, "Usage: %s [-a alignmentFile.paf] [-d RawLongReadsDir] [-k merSize] [-s minSupportForGoodRegions] [-l minLengthForGoodRegions] [-f freqThresholdForKMers] [-e maxError] [-p freqThresholdForKPers] [-m mode (0 for regions, 1 for cluster)] [-j threadsNb] \n\n", argv[0]);
				exit(EXIT_FAILURE);
        }
    }

    if (mode == 0) {
		goodRegionsGenK(alignmentFile, readsDir, merSize, minSupport, minRegionLength, thFreq, maxError, nbThreads, thPers);
	} else {
		clusterGenK(alignmentFile, readsDir, merSize, thFreq, maxError, nbThreads, thPers);
	}
}