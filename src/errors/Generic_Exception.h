/*
 * Generic_Exception.h
 *
 *  Created on: 7-set-2009
 *      Author: cdf
 */

#ifndef GENERIC_EXCEPTION_H_
#define GENERIC_EXCEPTION_H_

#include <exception>
#include <string>
using namespace std;

namespace errors {

class Generic_Exception: public std::exception {
public:
	Generic_Exception(const char* comment);
	virtual ~Generic_Exception() throw();
	virtual const char* what() const throw();
private:
	string _comment;
};


}

#endif /* GENERIC_EXCEPTION_H_ */
