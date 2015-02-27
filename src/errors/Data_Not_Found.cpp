#include "Data_Not_Found.h"

namespace errors
{

Data_Not_Found::~Data_Not_Found() throw()
{
	// Nothing to do
}

Data_Not_Found::Data_Not_Found(const char * comment) {
	stringstream ss;
	ss << "Data not found: " << comment << endl;
	_comment = ss.str();
}

const char* Data_Not_Found::what() const throw() {
	return _comment.c_str();
}

}
