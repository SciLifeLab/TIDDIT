#include "Data_Exception.h"

namespace errors {

Data_Exception::~Data_Exception() throw()
{
}

Data_Exception::Data_Exception(const string &comment) {
	_min = 0;
	_max = 0;
	_value = 0;
	_comment = string(comment);
	set_output();
}

Data_Exception::Data_Exception(const char *comment) {
	_min = 0;
	_max = 0;
	_value = 0;
	_comment = string(comment);
	set_output();
}

Data_Exception::Data_Exception(long int min, long int max, long int value) {
	_min = min;
	_max = max;
	_value = value;
	_comment = string("");
	set_output();
}

Data_Exception::Data_Exception(long int min, long int max, long int value, const string &comment) {
	_min = min;
	_max = max;
	_value = value;
	_comment = string(comment);
	set_output();
}

Data_Exception::Data_Exception(long int min, long int max, long int value, const char *comment) {
	_min = min;
	_max = max;
	_value = value;
	_comment = string(comment);
	set_output();
}

const char* Data_Exception::what() const throw() {
	return output.c_str();
}

void Data_Exception::add_comment(const string &str) {
	_comment.append(str);
	set_output();
}

void Data_Exception::add_comment(const char *str) {
	_comment.append(str);
	set_output();
}

void Data_Exception::set_output() {
	stringstream ss;
	ss << "Trying to access data at position \"" << _value
		<< "\" (range is: [" << _min << ".." << _max << "])";
	if (_comment.length() > 0) {
		ss << ": " << _comment;
	}
	output = ss.str();
}

}

