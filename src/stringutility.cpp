#include "stringutility.h"

void StringUtility::GetFirstArg(std::string &line, std::string &arg1)
{
	line = StringUtility::ReplaceString(line, "\t", " ");
	int istt = line.find_first_not_of(" ");
	if (istt == -1) 
	{
		arg1 = "";
	}
	else
	{
		line = line.substr(istt);
		int iend = line.find_first_of(" ");
		if (iend == -1) { iend = line.length(); }
		arg1 = line.substr(0,iend);
		line = line.substr(iend);
	}
}

std::string StringUtility::ReplaceString(std::string string1, std::string string2, std::string string3)
{
	std::string::size_type pos(string1.find(string2));
	while (pos != std::string::npos)
	{
		string1.replace(pos, string2.length(), string3);
		pos = string1.find(string2, pos+string3.length());
	}

	return string1;
}
