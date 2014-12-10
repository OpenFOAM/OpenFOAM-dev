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


#ifndef COCO_wmkdependPARSER_H__
#define COCO_wmkdependPARSER_H__

#include <iostream>
#include <string>
#include <list>
#include <set>

/*---------------------------------------------------------------------------*/



#include "wmkdependScanner.h"

namespace wmake {


/*---------------------------------------------------------------------------*\
                           Class Errors Declaration
\*---------------------------------------------------------------------------*/
//! Parser error handling
class Errors
{
public:
	int count;      //!< The number of errors detected

	//! Return a string describing the given error code.
	static std::wstring strerror(int n);

	Errors();               //!< Construct null - start with no errors
	virtual ~Errors();      //!< Destructor
	virtual void clear();   //!< Clear the error count

	//! Handle a general warning 'msg'
	virtual void Warning(const std::wstring& msg);
	//! Handle a general warning 'msg'
	virtual void Warning(int line, int col, const std::wstring& msg);
	//! Handle general error 'msg' (eg, a semantic error)
	virtual void Error(int line, int col, const std::wstring& msg);
	//! Handle syntax error 'n', uses strerror for the message, calls Error()
	virtual void SynErr(int line, int col, int n);
	//! Handle a general exception 'msg'
	virtual void Exception(const std::wstring& msg);

}; // Errors



/*---------------------------------------------------------------------------*\
                           Class Parser Declaration
\*---------------------------------------------------------------------------*/
//! A Coco/R Parser
class Parser
{
	enum {
		_EOF=0,
		_string=1,
		_sqstring=2,
		_package_name=3,
		_package_dir=4,
		maxT = 10    //<! max term (w/o pragmas)
	};
	static const int minErrDist = 2; //!< min. distance before reporting errors

	Token *dummyToken;
	bool deleteErrorsDestruct_; //!< delete the 'errors' member in destructor
	int  errDist;

	void SynErr(int n);         //!< Handle syntax error 'n'
	void Get();
	void Expect(int n);
	bool StartOf(int s);
	void ExpectWeak(int n, int follow);
	bool WeakSeparator(int n, int syFol, int repFol);

public:
	Scanner *scanner;
	Errors  *errors;

	Token *t;                   //!< last recognized token
	Token *la;                  //!< lookahead token

private:

    //! Set of (java) directories already visited
    static std::set<std::string> visitedDirs_;

    //! Replace all '.' with '/'
    static void dotToSlash(std::string& name);

    //! Import (java) directories
    static void importDir(const std::string& dirName);

    //! Import (java) file
    static void importFile(const std::string& name);

public:
    //! Set of files already visited
    static std::set<std::string> visitedFiles;

    //! Include directories to search
    static std::list<std::string> includeDirs;

    //! The name of the top-level source file
    static std::string sourceFile;

    //! The name of the top-level dep file
    static std::string depFile;

    //! Add directory to list of visited dirs, thus effectively ignoring it
    static void ignoreDir(const std::string& name);

    //! Include file
    static void includeFile(const std::string& name);

/*---------------------------------------------------------------------------*/

	//! Construct for the specified scanner
	/*!
	 * Use the default error handling, or optionally provide an error
	 * handler, which will not be deleted upon destruction.
	 */
	Parser(Scanner* scan, Errors* err = 0);
	~Parser();
	void Parse();                          //!< Execute the parse operation
	void SemErr(const std::wstring& msg);  //!< Handle semantic error
	bool isUTF8() const;   //!< Return true if scanner buffer is UTF8

	void wmkdepend();

}; // end Parser

} // End namespace

#endif // COCO_wmkdependPARSER_H__

// ************************************************************************* //
