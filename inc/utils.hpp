#ifndef _HG_UTILS_HPP_
#define _HG_UTILS_HPP_

#include <string>

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

	bool operator<(const alignment_t& a2) const;
};

struct consensus_t {
	std::string readId;
	int startPos;
	bool strand;
	int relPos;

	bool operator<(const consensus_t& c2) const;
};


namespace utils {

struct reverseComplement_t {
    
    reverseComplement_t();
   
    std::string operator()(std::string seq) const;

    private:
        char rev_comp_tab['T' + 1] = {0};
};

const static reverseComplement_t reverseComplement;
}

std::string getCannonicalSeq(std::string seq);

#endif //_HG_UTILS_HPP_
