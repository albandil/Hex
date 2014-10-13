//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
//                                                                                   //
//                       / /   / /    __    \ \  / /                                 //
//                      / /__ / /   / _ \    \ \/ /                                  //
//                     /  ___  /   | |/_/    / /\ \                                  //
//                    / /   / /    \_\      / /  \ \                                 //
//                                                                                   //
//                                                                                   //
//  Copyright (c) 2014, Jakub Benda, Charles University in Prague                    //
//                                                                                   //
// MIT License:                                                                      //
//                                                                                   //
//  Permission is hereby granted, free of charge, to any person obtaining a          //
// copy of this software and associated documentation files (the "Software"),        //
// to deal in the Software without restriction, including without limitation         //
// the rights to use, copy, modify, merge, publish, distribute, sublicense,          //
// and/or sell copies of the Software, and to permit persons to whom the             //
// Software is furnished to do so, subject to the following conditions:              //
//                                                                                   //
//  The above copyright notice and this permission notice shall be included          //
// in all copies or substantial portions of the Software.                            //
//                                                                                   //
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS          //
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE       //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, //
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF         //
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.  //
//                                                                                   //
//  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //

#include <iostream>
#include <vector>
#include <string>

#include "../../src/arrays.h"
#include "../../src/arrays.cpp"
#include "../../src/hdffile.cpp"

int main(int argc, char *argv[])
{
	H5::Exception::dontPrint();
	
	std::vector<std::string> files;
	bool compress = false;
	
	for (int iarg = 1; iarg < argc; iarg++)
	{
		if (strcmp(argv[iarg],"--compress") == 0)
			compress = true;
		else
			files.push_back(std::string(argv[iarg]));
	}
	
	if (files.size() != 2)
	{
		std::cout << "\nUsage:\n\thdf2hdf [--compress] <file1> <file2>\n";
		return -1;
	}
	
	rArray a;
	
	if (not a.hdfload(files[0].c_str()))
	{
		std::cout << "Can't open the file \"" << argv[1] << "\"\n";
		return -1;
	}
	
	if (not a.hdfsave(files[1].c_str(), compress))
	{
		std::cout << "Error while writing file \"" << argv[2] << "\"\n";
		return -1;
	}
	
	return 0;
}
