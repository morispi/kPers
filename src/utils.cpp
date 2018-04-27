#include "utils.hpp"

bool alignment_t::operator<(const alignment_t& a2) const {
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

bool consensus_t::operator<(const consensus_t& c2) const {
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

utils::reverseComplement_t::reverseComplement_t()
{
    rev_comp_tab['A'] = 'T';
    rev_comp_tab['T'] = 'A';
    rev_comp_tab['C'] = 'G';
    rev_comp_tab['G'] = 'C';
    rev_comp_tab['a'] = 't';
    rev_comp_tab['t'] = 'a';
    rev_comp_tab['c'] = 'g';
    rev_comp_tab['g'] = 'c';
}


std::string utils::reverseComplement_t::operator()(std::string seq) const  {

    auto first = seq.begin(), last = seq.end();

    while(true)
	if(first == last || first == --last)
	{
	    if(seq.length() % 2)
		*first = this->rev_comp_tab[*first];
	    return seq;
	}
	else
	{
	    *first = this->rev_comp_tab[*first];
	    *last = this->rev_comp_tab[*last];
	    std::iter_swap(first, last);
	    ++first;
	}

    return seq;
}	

std::string getCannonicalSeq(std::string seq) {
	std::string revComp = utils::reverseComplement(seq);
	return seq < revComp ? seq : revComp;
}

