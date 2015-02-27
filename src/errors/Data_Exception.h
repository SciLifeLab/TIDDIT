#ifndef DATA_EXCEPTION_H_
#define DATA_EXCEPTION_H_

#include <exception>
#include <sstream>
using namespace std;

namespace errors {

class Data_Exception: public std::exception
{
public:
	virtual ~Data_Exception() throw();
	Data_Exception(const string &comment);
	Data_Exception(const char *comment);
	Data_Exception(long int min, long int max, long int value);
	Data_Exception(long int min, long int max, long int value, const string &comment);
	Data_Exception(long int min, long int max, long int value, const char *comment);
	virtual const char* what() const throw();
	void add_comment(const string &str);
	void add_comment(const char *str);
private:
	long int _min;
	long int _max;
	long int _value;
	string _comment;
	string output;

	void set_output();
};

}

#endif /*DATA_EXCEPTION_H_*/
