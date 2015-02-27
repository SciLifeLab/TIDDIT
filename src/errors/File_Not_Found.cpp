#include "File_Not_Found.h"

namespace errors
{

File_Not_Found::~File_Not_Found() throw()
{
}

File_Not_Found::File_Not_Found(const char* comment)
{
	stringstream ss;
	ss << "File not found: " << comment << endl;
	_comment = ss.str();
}

File_Not_Found::File_Not_Found(const string &comment) {
	stringstream ss;
	ss << "File not found: " << comment << endl;
	_comment = ss.str();
}

File_Not_Found::File_Not_Found(const string* comment) {
	stringstream ss;
	ss << "File not found: " << *comment << endl;
	_comment = ss.str();
}

const char* File_Not_Found::what() const throw() {
	return _comment.c_str();
}


}
