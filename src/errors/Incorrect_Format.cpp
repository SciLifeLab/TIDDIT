#include "Incorrect_Format.h"

namespace errors
{

Incorrect_Format::~Incorrect_Format() throw()
{
}

Incorrect_Format::Incorrect_Format(const string &comment) {
	stringstream ss;
	ss << "Incorrect format: " << comment << endl;
	_comment = ss.str();
}

Incorrect_Format::Incorrect_Format(const string* comment) {
	stringstream ss;
	ss << "Incorrect format: " << *comment << endl;
	_comment = ss.str();
}

Incorrect_Format::Incorrect_Format(const char* comment) {
	stringstream ss;
	ss << "Incorrect format: " << comment << endl;
	_comment = ss.str();
}

const char* Incorrect_Format::what() const throw() {
	return _comment.c_str();
}

}
