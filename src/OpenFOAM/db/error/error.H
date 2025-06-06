/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

Class
    Foam::error

Description
    Class to handle errors and exceptions in a simple, consistent stream-based
    manner.

    The error class is globally instantiated with a title string. Errors,
    messages and other data are piped to the messageStream class in the
    standard manner.  Manipulators are supplied for exit and abort which may
    terminate the program or throw an exception depending on whether the
    exception handling has been switched on (off by default).

Usage
    \code
        error << "message1" << "message2" << FoamDataType << exit(errNo);
        error << "message1" << "message2" << FoamDataType << abort();
    \endcode

SourceFiles
    error.C

\*---------------------------------------------------------------------------*/

#ifndef error_H
#define error_H

#include "OStringStream.H"
#include "messageStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class error;
Ostream& operator<<(Ostream&, const error&);


/*---------------------------------------------------------------------------*\
                            Class error Declaration
\*---------------------------------------------------------------------------*/

class error
:
    public std::exception,
    public messageStream
{

protected:

    // Protected data

        string functionName_;
        string sourceFileName_;
        label sourceFileLineNumber_;

        bool abort_;

        bool throwExceptions_;

        OStringStream messageStream_;


public:

    // Constructors

        //- Construct from title string
        error(const string& title);


    // Member Functions

        string message() const;

        const string& functionName() const
        {
            return functionName_;
        }

        const string& sourceFileName() const
        {
            return sourceFileName_;
        }

        label sourceFileLineNumber() const
        {
            return sourceFileLineNumber_;
        }

        void throwExceptions()
        {
            throwExceptions_ = true;
        }

        void dontThrowExceptions()
        {
            throwExceptions_ = false;
        }

        //- Convert to OSstream
        //  Prints basic message and returns OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Convert to OSstream
        //  Prints basic message and returns OSstream for further info.
        OSstream& operator()
        (
            const string& functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber = 0
        );

        //- Explicitly convert to OSstream for << operations
        OSstream& operator()();

        //- Convert to OSstream for << operations
        operator OSstream&()
        {
            return this->operator()();
        }

        //- Create and return a dictionary
        operator dictionary() const;

        //- Helper function to print a stack (if OpenFOAM IO not yet
        //  initialised)
        static void safePrintStack(std::ostream&);

        //- Helper function to print a stack
        static void printStack(Ostream&);

        //- Exit : can be called for any error to exit program.
        //  Prints stack before exiting.
        void exit(const int errNo = 1);

        //- Abort : used to stop code for fatal errors.
        //  Prints stack before exiting.
        void abort();


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const error&);
};


// Forward declaration of friend functions and operators
class IOerror;
Ostream& operator<<(Ostream&, const IOerror&);


/*---------------------------------------------------------------------------*\
                       Class IOerrorLocation Declaration
\*---------------------------------------------------------------------------*/

class IOerrorLocation
{
    // Private Data

        //- File name
        string ioFileName_;

        //- Start line number
        label ioStartLineNumber_;

        //- End line number
        label ioEndLineNumber_;

        //- Is file global
        bool ioGlobal_;


public:

    // Constructors

        //- Construct null
        IOerrorLocation();

        //- Construct from components
        IOerrorLocation
        (
            const string& ioFileName,
            const label ioStartLineNumber = -1,
            const label ioEndLineNumber = -1,
            const bool ioGlobal = false
        );

        //- Construct from a stream
        IOerrorLocation(const IOstream&);

        //- Construct from a dictionary
        IOerrorLocation(const dictionary&);


    // Member Functions

        //- Access the file name
        inline const string& ioFileName() const
        {
            return ioFileName_;
        }

        //- Access the start line number
        inline label ioStartLineNumber() const
        {
            return ioStartLineNumber_;
        }

        //- Access the end line number
        inline label ioEndLineNumber() const
        {
            return ioEndLineNumber_;
        }

        //- Assess to the global/local file switch
        inline bool ioGlobal() const
        {
            return ioGlobal_;
        }
};


/*---------------------------------------------------------------------------*\
                           Class IOerror Declaration
\*---------------------------------------------------------------------------*/

class IOerror
:
    public error,
    public IOerrorLocation
{
public:

    // Constructors

        //- Construct from title string
        IOerror(const string& title);


    // Member Functions

        //- Convert to OSstream
        //  Prints basic message and returns OSstream for further info.
        OSstream& operator()
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const IOerrorLocation& location
        );

        //- Print basic message and exit. Uses cerr if streams not constructed
        //  yet (at startup). Use in startup parsing instead of FatalError.
        static void SafeFatalIOError
        (
            const char* functionName,
            const char* sourceFileName,
            const int sourceFileLineNumber,
            const IOstream&,
            const string& msg
        );

        //- Create and return a dictionary
        operator dictionary() const;

        //- Exit : can be called for any error to exit program
        void exit(const int errNo = 1);

        //- Abort : used to stop code for fatal errors
        void abort();


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const IOerror&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error declarations: defined in error.C

extern error   FatalError;
extern IOerror FatalIOError;

// Template argument dependent "false" for static_assert in templated functions
template <class... T>
constexpr bool False = false;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Convenience macros to add the file name and line number to the function name

//- Report an error message using Foam::FatalError
//  for functionName in file __FILE__ at line __LINE__
#define FatalErrorIn(functionName)                                             \
    ::Foam::FatalError((functionName), __FILE__, __LINE__)

//- Report an error message using Foam::FatalError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
#define FatalErrorInFunction FatalErrorIn(FUNCTION_NAME)


//- Report an error message using Foam::FatalIOError
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define FatalIOErrorIn(functionName, ios)                                      \
    ::Foam::FatalIOError((functionName), __FILE__, __LINE__, (ios))

//- Report an error message using Foam::FatalIOError
//  for FUNCTION_NAME in file __FILE__ at line __LINE__
//  for a particular IOstream
#define FatalIOErrorInFunction(ios) FatalIOErrorIn(FUNCTION_NAME, ios)


//- Report an error message using Foam::FatalIOError
//  (or cerr if FatalIOError not yet constructed)
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SafeFatalIOErrorIn(functionName, ios, msg)                             \
    ::Foam::IOerror::SafeFatalIOError                                          \
    ((functionName), __FILE__, __LINE__, (ios), (msg))

//- Report an error message using Foam::FatalIOError
//  (or cerr if FatalIOError not yet constructed)
//  for functionName in file __FILE__ at line __LINE__
//  for a particular IOstream
#define SafeFatalIOErrorInFunction(ios, msg)                                   \
    SafeFatalIOErrorIn(FUNCTION_NAME, ios, msg)


//- Issue a FatalErrorIn for a function not currently implemented.
//  The functionName is printed and then abort is called.
//
//  This macro can be particularly useful when methods must be defined to
//  complete the interface of a derived class even if they should never be
//  called for this derived class.
#define notImplemented(functionName)                                           \
    FatalErrorIn(functionName)                                                 \
        << "Not implemented" << ::Foam::abort(::Foam::FatalError);

//- Issue a FatalErrorIn for a function not currently implemented.
//  The FUNCTION_NAME is printed and then abort is called.
//
//  This macro can be particularly useful when methods must be defined to
//  complete the interface of a derived class even if they should never be
//  called for this derived class.
#define NotImplemented notImplemented(FUNCTION_NAME)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "errorManip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
