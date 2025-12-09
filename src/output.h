#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <iostream>
#include <string>
#include <cstdio>

/* class for output on screen */
class Output
{
	public:
		static void ShowError(const char*);
		static void ShowWarning(const char*);
		static void ShowMessage(const char*);
};

#endif // _OUTPUT_H
