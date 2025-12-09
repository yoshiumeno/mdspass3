#ifndef STRINGUTILITY_H
#define STRINGUTILITY_H


#include <iostream>
#include <string>
#include <sstream>

class StringUtility
{
	public:

	static void GetFirstArg(std::string&, std::string&);
	static std::string ReplaceString(std::string string1, std::string string2, std::string string3);
};

#endif // STRINGUTILITY_H
