#ifndef TYPES_H_
#define TYPES_H_

#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdlib>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

#ifdef INLINE_DISABLED
#define INLINE
#else
#define INLINE inline
#endif



#include <boost/program_options.hpp>
namespace po = boost::program_options;


#include <string>
using namespace std;


#define DEFAULT_CHANNEL std::cout
#define ERROR_CHANNEL std::cerr

#define VERBOSE_CHANNEL std::cerr
#define DEBUG_CHANNEL std::cerr
#define DEFAULT_CHANNEL std::cout

static inline std::string package_description() {
	std::string line(PACKAGE_NAME);
	line.append(" version ");
	line.append(PACKAGE_VERSION);
	return line;
}

#define text_delimitator '$'
#define masked_base 'N'
#define match_character '.'
#define removed_character '-'

typedef std::size_t t_size; ///< Type for Vectors or similar classes

//#define LONG_LENGTH

#ifdef LONG_LENGTH
	typedef long int t_length; ///< Type for TEXT
	#define MAX_T_LENGTH LONG_MAX ///< Max value for TEXT
	#define SA_NOT_FOUND -1
#else
	typedef int t_length; ///< Type for TEXT
	#define MAX_T_LENGTH INT_MAX ///< Max value for TEXT
	#define SA_NOT_FOUND -1
#endif

typedef unsigned short int t_char; ///< Type for char conversion
typedef unsigned short int t_errors; ///< Type for ERRORS
typedef t_errors t_errors_delta ;
typedef unsigned int t_pattern_length; ///< Type for PATTERN
typedef double t_quality; ///< Type for QUALITY VALUES

typedef std::vector<short unsigned int> t_quality_vector; ///< Vector with quality values

typedef std::string t_edit_string; ///< For future expansion

enum t_strand { forward_strand, reverse_strand, unknown_strand }; ///< Just an enumeration of possible strands

enum t_alignment { unique_alignment, multiple_alignments, alignments_not_found, quality_discarded, low_complexity, contamination, unknown_alignment };

enum t_masked { masked, not_masked };

enum Feature {LOW_COVERAGE_AREA, HIGH_COVERAGE_AREA, LOW_NORMAL_AREA, HIGH_NORMAL_AREA, HIGH_SINGLE_AREA, HIGH_SPANNING_AREA, HIGH_OUTIE_AREA, COMPRESSION_AREA, STRECH_AREA, TOTAL};



/**
 * Conversion from strand type to char type
 */
static inline char strand_to_char(const t_strand s) {
	if (s == unknown_strand)
		return '.';
	else
		return (s == forward_strand) ? '+' : '-';
}

/**
 * Conversion from char type to strand type
 */
static inline t_strand char_to_strand(const char c) {
	if (c == '+' or c == 'F')
		return forward_strand;
	else if (c == '-' or c == 'R')
		return reverse_strand;
	else
		return unknown_strand;
}

#define unique_alignment_char 'U'
#define multiple_alignments_char 'M'
#define alignments_not_found_char 'N'
#define quality_discarded_char 'Q'
#define low_complexity_char 'L'
#define contamination_char 'C'
#define unknown_alignment_to_char '?'

#define unique_alignment_string "U"
#define multiple_alignments_string "M"
#define alignments_not_found_string "NF"
#define quality_discarded_string "QD"
#define low_complexity_string "LC"
#define contamination_string "C"
#define unknown_alignment_to_string "?"

/**
 * Conversion from alignment type to char type
 */
static inline char alignment_to_char(const t_alignment a) {
	switch (a) {
	case unique_alignment :		return unique_alignment_char;
	case multiple_alignments :	return multiple_alignments_char;
	case alignments_not_found :	return alignments_not_found_char;
	case quality_discarded :	return quality_discarded_char;
	case low_complexity :		return low_complexity_char;
	case contamination :		return contamination_char;
	case unknown_alignment :	return unknown_alignment_to_char;
	}
	return unknown_alignment_to_char;
}

/**
 * Conversion from char type to alignment type
 */
static inline t_alignment char_to_alignment(char c) {
	switch (c) {
	case unique_alignment_char :		return unique_alignment;
	case 'R': // backward compatibility
	case multiple_alignments_char :		return multiple_alignments;
	case alignments_not_found_char :	return alignments_not_found;
	case quality_discarded_char :		return quality_discarded;
	case low_complexity_char :			return low_complexity;
	case contamination_char : 			return contamination;
	default : return unknown_alignment;
	}
}

/**
 * Conversion from alignment type to string type
 */
static inline std::string alignment_to_string(const t_alignment a) {
	switch (a) {
	case unique_alignment :		return std::string(unique_alignment_string);
	case multiple_alignments :	return std::string(multiple_alignments_string);
	case alignments_not_found :	return std::string(alignments_not_found_string);
	case quality_discarded :	return std::string(quality_discarded_string);
	case low_complexity :		return std::string(low_complexity_string);
	case contamination :		return std::string(contamination_string);
	case unknown_alignment :	return std::string(unknown_alignment_to_string);
	}
	return std::string(unknown_alignment_to_string);
}

static inline char reverse_complement_standalone_char(const char c) {
	switch (c) {
	case 'A' : return 'T';
	case 'T' : return 'A';
	case 'C' : return 'G';
	case 'G' : return 'C';
	case 'U' : return 'A';
	case 'R' : return 'Y';
	case 'Y' : return 'R';
	case 'M' : return 'K';
	case 'K' : return 'M';
	case 'W' : return 'S';
	case 'S' : return 'W';
	case 'B' : return 'V';
	case 'V' : return 'B';
	case 'D' : return 'H';
	case 'H' : return 'D';
	//case 'N' : return 'N';
	//case 'X' : return 'X';
	default : return c;
	}
}

static inline std::string reverse_complement_standalone_str_length(const char *str, size_t length) {
	std::string reverse;
	t_size i = length;
	while (i > 0)
		reverse.push_back(reverse_complement_standalone_char(str[--i]));
	return reverse;
}

static inline std::string reverse_complement_standalone_str(const char *str) {
	return reverse_complement_standalone_str_length(str,strlen(str));
}

static inline std::string reverse_complement_standalone_str(const std::string & str) {
	return reverse_complement_standalone_str_length(str.c_str(),str.length());
}


static inline std::string reverse_standalone_str_length(const char *str, size_t length) {
	std::string reverse;
	t_size i = length;
	while (i > 0)
		reverse.push_back(str[--i]);
	return reverse;
}

static inline std::string reverse_standalone_str(const char *str) {
	return reverse_standalone_str_length(str,strlen(str));
}

static inline std::string tolower(const std::string & old_string) {
	std::string lower_string(old_string);
	for (std::string::iterator iter = lower_string.begin(); iter != lower_string.end(); iter++) {
		*iter = tolower(*iter);
	}
	return lower_string;
}

static inline std::string toupper(const std::string & old_string){
	std::string upper_string(old_string);
	for (std::string::iterator iter = upper_string.begin(); iter != upper_string.end(); iter++) {
		*iter = toupper(*iter);
	}
	return upper_string;
}

#endif /*TYPES_H_*/
