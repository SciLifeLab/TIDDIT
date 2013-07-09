#ifndef FILE_NOT_FOUND_H_
#define FILE_NOT_FOUND_H_

#include <exception>
#include <sstream>
using namespace std;

namespace errors
{

class File_Not_Found: public std::exception
{
public:
	virtual ~File_Not_Found() throw();
	File_Not_Found(const char* comment);
	File_Not_Found(const string &comment);
	File_Not_Found(const string* comment);
	virtual const char* what() const throw();

private:
	string _comment;
};

}

#endif /*FILE_NOT_FOUND_H_*/
