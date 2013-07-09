#ifndef DATA_NOT_FOUND_H_
#define DATA_NOT_FOUND_H_

#include <exception>
#include <sstream>
using namespace std;

namespace errors
{

class Data_Not_Found: public std::exception
{
public:
	virtual ~Data_Not_Found() throw();
	Data_Not_Found(const char* comment);
	virtual const char* what() const throw();
private:
	string _comment;
};

}

#endif /*DATA_NOT_FOUND_H_*/
