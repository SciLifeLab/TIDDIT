#ifndef INCORRECT_FORMAT_H_
#define INCORRECT_FORMAT_H_

#include <exception>
#include <sstream>
using namespace std;

namespace errors
{

class Incorrect_Format: public std::exception
{
public:
	virtual ~Incorrect_Format() throw();
	Incorrect_Format(const string &comment);
	Incorrect_Format(const string* comment);
	Incorrect_Format(const char* comment);
	virtual const char* what() const throw();
private:
	string _comment;
};

}

#endif /*INCORRECT_FORMAT_H_*/
