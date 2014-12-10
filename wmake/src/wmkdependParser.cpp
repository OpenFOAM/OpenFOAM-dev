/*---------------------------------*- C++ -*---------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

@file wmkdependParser.atg

Description
    An attributed Coco/R grammar to parse C/C++, Fortran and Java files
    for include and import statements.

SourceFiles
    generated

\*---------------------------------------------------------------------------*/
// This file was generated with Coco/R C++ (10 Mar 2010)
// http://www.ssw.uni-linz.ac.at/coco/
// with these defines:
//     - FORCE_UTF8


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cwchar>
#include <sstream>

#include "wmkdependParser.h"

namespace wmake {


#include <sys/types.h>
#include <dirent.h>

std::set<std::string> Parser::visitedDirs_;

std::set<std::string> Parser::visitedFiles;
std::list<std::string> Parser::includeDirs;
std::string Parser::sourceFile;
std::string Parser::depFile;


void Parser::dotToSlash(std::string& name)
{
    std::string::size_type start = 0;

    while ((start = name.find('.', start)) != std::string::npos)
    {
        name.replace(start, 1, 1, '/');
        start++;
    }
}


void Parser::ignoreDir(const std::string& name)
{
    visitedDirs_.insert(name);
}


void Parser::includeFile(const std::string& name)
{
    if (!visitedFiles.insert(name).second)
    {
        return;   // already existed (did not insert)
    }

    // use stdio and buffering within Coco/R -- (faster)
    FILE *fh = fopen(name.c_str(), "r");
    if (fh)
    {
        std::cout << depFile << ": " << name << "\n";
    }
    else
    {
        for
        (
            std::list<std::string>::const_iterator iter = includeDirs.begin();
            iter != includeDirs.end();
            ++iter
        )
        {
            const std::string pathName = *iter + name;

            fh = fopen(pathName.c_str(), "r");
            if (fh)
            {
                std::cout << depFile << ": " << pathName << "\n";
                break;
            }
        }
    }

    if (fh)
    {
        Scanner scanner(fh);
        Parser  parser(&scanner);

        parser.Parse();
        fclose(fh);
    }
    else
    {
        fwprintf
        (
            stderr,
            L"could not open file %s for source file %s\n",
            name.c_str(), sourceFile.c_str()
        );

        // only report the first occurance
        visitedFiles.insert(name);
    }
}


void Parser::importFile(const std::string& name)
{
    // check if a globbed form was already visited
    std::string::size_type dotPos = name.find('.');
    if (dotPos != std::string::npos)
    {
        std::string dirGlob = name.substr(0, dotPos);
        dirGlob += ".*";

        if (visitedDirs_.find(dirGlob) != visitedDirs_.end())
        {
            return;
        }
    }

    std::string javaFileName = name;

    dotToSlash(javaFileName);
    javaFileName += ".java";

    includeFile(javaFileName);
}


void Parser::importDir(const std::string& name)
{
    if (!visitedDirs_.insert(name).second)
    {
        return;   // already existed (did not insert)
    }

    std::string dirName = name;
    dotToSlash(dirName);

    DIR *source = opendir(dirName.c_str());

    if (source)
    {
        struct dirent *list;

        // Read and parse all the entries in the directory
        while ((list = readdir(source)) != NULL)
        {
            const char* ext = strstr(list->d_name, ".java");

            // avoid matching on something like '.java~'
            if (ext && strlen(ext) == 5)
            {
                std::string pathName = dirName + list->d_name;
                includeFile(pathName);
            }
        }

        closedir(source);
    }
    else
    {
        fwprintf
        (
            stderr,
            L"could not open directory %s\n",
            dirName.c_str()
        );
        return;
    }
}




//! @cond fileScope
//
//  Create by copying str - only used locally
inline static wchar_t* coco_string_create(const wchar_t* str)
{
	const int len = wcslen(str);
	wchar_t* dst = new wchar_t[len + 1];
	wcsncpy(dst, str, len);
	dst[len] = 0;
	return dst;
}


// Free storage and nullify the argument
inline static void coco_string_delete(wchar_t* &str)
{
	delete[] str;
	str = NULL;
}
//
//! @endcond


// ----------------------------------------------------------------------------
// Parser Implementation
// ----------------------------------------------------------------------------

void Parser::SynErr(int n)
{
	if (errDist >= minErrDist) errors->SynErr(la->line, la->col, n);
	errDist = 0;
}


void Parser::SemErr(const std::wstring& msg)
{
	if (errDist >= minErrDist) errors->Error(t->line, t->col, msg);
	errDist = 0;
}


bool Parser::isUTF8() const
{
	return scanner && scanner->buffer && scanner->buffer->isUTF8();
}


void Parser::Get()
{
	for (;;)
	{
		t = la;
		la = scanner->Scan();
		if (la->kind <= maxT)
		{
			++errDist;
			break;
		}
		if (dummyToken != t)
		{
			dummyToken->kind = t->kind;
			dummyToken->pos = t->pos;
			dummyToken->col = t->col;
			dummyToken->line = t->line;
			dummyToken->next = NULL;
			coco_string_delete(dummyToken->val);
			dummyToken->val = coco_string_create(t->val);
			t = dummyToken;
		}
		la = t;
	}
}


void Parser::Expect(int n)
{
	if (la->kind == n)
	{
		Get();
	}
	else
	{
		SynErr(n);
	}
}


void Parser::ExpectWeak(int n, int follow)
{
	if (la->kind == n)
	{
		Get();
	}
	else
	{
		SynErr(n);
		while (!StartOf(follow))
		{
			Get();
		}
	}
}


bool Parser::WeakSeparator(int n, int syFol, int repFol)
{
	if (la->kind == n)
	{
		Get();
		return true;
	}
	else if (StartOf(repFol))
	{
		return false;
	}
	else
	{
		SynErr(n);
		while (!(StartOf(syFol) || StartOf(repFol) || StartOf(0)))
		{
			Get();
		}
		return StartOf(syFol);
	}
}


void Parser::wmkdepend()
{
	while (StartOf(1)) {
		if (la->kind == 5) {
			Get();
			if (la->kind == 6) {
				Get();
				if (la->kind == 1) {
					Get();
					if (isUTF8())
					{
					    includeFile(t->toStringUTF8(1, t->length()-2));
					}
					else
					{
					    includeFile(t->toString(1, t->length()-2));
					}

				}
			}
			if (StartOf(2)) {
				Get();
				while (StartOf(3)) {
					Get();
				}
			}
			Expect(7);
		} else if (la->kind == 6) {
			Get();
			if (la->kind == 2) {
				Get();
				if (isUTF8())
				{
				    includeFile(t->toStringUTF8(1, t->length()-2));
				}
				else
				{
				    includeFile(t->toString(1, t->length()-2));
				}

			}
			if (StartOf(4)) {
				Get();
				while (StartOf(3)) {
					Get();
				}
			}
			Expect(7);
		} else if (la->kind == 8) {
			Get();
			if (la->kind == 4) {
				Get();
				if (isUTF8())
				{
				    importDir(t->toStringUTF8());
				}
				else
				{
				    importDir(t->toString());
				}

			} else if (la->kind == 3) {
				Get();
				if (isUTF8())
				{
				    importFile(t->toStringUTF8());
				}
				else
				{
				    importFile(t->toString());
				}

			} else SynErr(11);
			Expect(9);
			if (StartOf(3)) {
				Get();
				while (StartOf(3)) {
					Get();
				}
			}
			Expect(7);
		} else {
			if (StartOf(5)) {
				Get();
				while (StartOf(3)) {
					Get();
				}
			}
			Expect(7);
		}
	}
}



void Parser::Parse()
{
	t = NULL;
	// might call Parse() twice
	if (dummyToken) {
		coco_string_delete(dummyToken->val);
		delete dummyToken;
	}
	dummyToken = new Token(coco_string_create(L"Dummy Token"));
	la = dummyToken;
	Get();
	wmkdepend();
	Expect(0);  // expect end-of-file automatically added
}


Parser::Parser(Scanner* scan, Errors* err)
:
	dummyToken(NULL),
	deleteErrorsDestruct_(!err),
	errDist(minErrDist),
	scanner(scan),
	errors(err),
	t(NULL),
	la(NULL)
{
	if (!errors)    // add in default error handling
	{
		errors = new Errors();
	}
	// user-defined initializations:
}


bool Parser::StartOf(int s)
{
	const bool T = true;
	const bool x = false;

	static const bool set[6][12] =
	{
		{T,x,x,x, x,x,x,x, x,x,x,x},
		{x,T,T,T, T,T,T,T, T,T,T,x},
		{x,x,T,T, T,T,x,x, T,T,T,x},
		{x,T,T,T, T,T,T,x, T,T,T,x},
		{x,T,x,T, T,T,T,x, T,T,T,x},
		{x,T,T,T, T,x,x,x, x,T,T,x}
	};

	return set[s][la->kind];
}


Parser::~Parser()
{
	if (deleteErrorsDestruct_) { delete errors; } // delete default error handling
	if (dummyToken) {
		coco_string_delete(dummyToken->val);
		delete dummyToken;
	}
	// user-defined destruction:
}


// ----------------------------------------------------------------------------
// Errors Implementation
// ----------------------------------------------------------------------------

Errors::Errors()
:
	count(0)
{}


Errors::~Errors()
{}


void Errors::clear()
{
	count = 0;
}


std::wstring Errors::strerror(int n)
{
	switch (n) {
		case 0: return L"EOF expected"; break;
		case 1: return L"string expected"; break;
		case 2: return L"sqstring expected"; break;
		case 3: return L"package_name expected"; break;
		case 4: return L"package_dir expected"; break;
		case 5: return L"\"#\" expected"; break;
		case 6: return L"\"include\" expected"; break;
		case 7: return L"\"\\n\" expected"; break;
		case 8: return L"\"import\" expected"; break;
		case 9: return L"\";\" expected"; break;
		case 10: return L"??? expected"; break;
		case 11: return L"invalid wmkdepend"; break;
		default:
		{
			// std::wostringstream buf;  (this typedef might be missing)
			std::basic_ostringstream<wchar_t> buf;
			buf << "error " << n;
			return buf.str();
		}
		break;
	}
}


void Errors::Warning(const std::wstring& msg)
{
	fwprintf(stderr, L"%ls\n", msg.c_str());
}


void Errors::Warning(int line, int col, const std::wstring& msg)
{
	fwprintf(stderr, L"-- line %d col %d: %ls\n", line, col, msg.c_str());
}


void Errors::Error(int line, int col, const std::wstring& msg)
{
	fwprintf(stderr, L"-- line %d col %d: %ls\n", line, col, msg.c_str());
	count++;
}


void Errors::SynErr(int line, int col, int n)
{
	this->Error(line, col, this->strerror(n));
}


void Errors::Exception(const std::wstring& msg)
{
	fwprintf(stderr, L"%ls", msg.c_str());
	::exit(1);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace

// ************************************************************************* //
