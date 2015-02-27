/*
 * Generic_Exception.cpp
 *
 *  Created on: 7-set-2009
 *      Author: cdf
 */

#include "Generic_Exception.h"

namespace errors {

Generic_Exception::Generic_Exception(const char* comment) {
	_comment = string(comment);
}

Generic_Exception::~Generic_Exception() throw() {
	// Nothing to do
}

const char* Generic_Exception::what() const throw() {
	return _comment.c_str();
}


}
