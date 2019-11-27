/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "OStringStream.H"
#include "fileName.H"
#include "dictionary.H"
#include "jobInfo.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::IOerror::IOerror(const string& title)
:
    error(title),
    ioFileName_("unknown"),
    ioStartLineNumber_(-1),
    ioEndLineNumber_(-1)
{}


Foam::IOerror::IOerror(const dictionary& errDict)
:
    error(errDict),
    ioFileName_(errDict.lookup("ioFileName")),
    ioStartLineNumber_(errDict.lookup<label>("ioStartLineNumber")),
    ioEndLineNumber_(errDict.lookup<label>("ioEndLineNumber"))
{}


Foam::IOerror::~IOerror() throw()
{}


Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const string& ioFileName,
    const label ioStartLineNumber,
    const label ioEndLineNumber
)
{
    error::operator()(functionName, sourceFileName, sourceFileLineNumber);
    ioFileName_ = ioFileName;
    ioStartLineNumber_ = ioStartLineNumber;
    ioEndLineNumber_ = ioEndLineNumber;

    return operator OSstream&();
}


Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const IOstream& ioStream
)
{
    return operator()
    (
        functionName,
        sourceFileName,
        sourceFileLineNumber,
        ioStream.name(),
        ioStream.lineNumber(),
        -1
    );
}


Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const dictionary& dict
)
{
    return operator()
    (
        functionName,
        sourceFileName,
        sourceFileLineNumber,
        dict.name(),
        dict.startLineNumber(),
        dict.endLineNumber()
    );
}


void Foam::IOerror::SafeFatalIOError
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const IOstream& ioStream,
    const string& msg
)
{
    if (jobInfo::constructed)
    {
        FatalIOError
        (
            functionName,
            sourceFileName,
            sourceFileLineNumber,
            ioStream
        )   << msg << Foam::exit(FatalIOError);
    }
    else
    {
        std::cerr
            << std::endl
            << "--> FOAM FATAL IO ERROR:" << std::endl
            << msg
            << std::endl
            << "file: " << ioStream.name()
            << " at line " << ioStream.lineNumber() << '.'
            << std::endl << std::endl
            << "    From function " << functionName
            << std::endl
            << "    in file " << sourceFileName
            << " at line " << sourceFileLineNumber << '.'
            << std::endl;
        ::exit(1);
    }
}


Foam::IOerror::operator Foam::dictionary() const
{
    dictionary errDict(error::operator dictionary());

    errDict.remove("type");
    errDict.add("type", word("Foam::IOerror"));

    errDict.add("ioFileName", ioFileName());
    errDict.add("ioStartLineNumber", ioStartLineNumber());
    errDict.add("ioEndLineNumber", ioEndLineNumber());

    return errDict;
}


void Foam::IOerror::exit(const int)
{
    if (!throwExceptions_ && jobInfo::constructed)
    {
        jobInfo_.add("FatalIOError", operator dictionary());
        jobInfo_.exit();
    }

    if (abort_)
    {
        abort();
    }

    if (Pstream::parRun())
    {
        Perr<< endl << *this << endl
            << "\nFOAM parallel run exiting\n" << endl;
        Pstream::exit(1);
    }
    else
    {
        if (throwExceptions_)
        {
            // Make a copy of the error to throw
            IOerror errorException(*this);

            // Rewind the message buffer for the next error message
            messageStreamPtr_->rewind();

            throw errorException;
        }
        else
        {
            Perr<< endl << *this << endl
                << "\nFOAM exiting\n" << endl;
            ::exit(1);
        }
    }
}


void Foam::IOerror::abort()
{
    if (!throwExceptions_ && jobInfo::constructed)
    {
        jobInfo_.add("FatalIOError", operator dictionary());
        jobInfo_.abort();
    }

    if (abort_)
    {
        Perr<< endl << *this << endl
            << "\nFOAM aborting (FOAM_ABORT set)\n" << endl;
        printStack(Perr);
        ::abort();
    }

    if (Pstream::parRun())
    {
        Perr<< endl << *this << endl
            << "\nFOAM parallel run aborting\n" << endl;
        printStack(Perr);
        Pstream::abort();
    }
    else
    {
        if (throwExceptions_)
        {
            // Make a copy of the error to throw
            IOerror errorException(*this);

            // Rewind the message buffer for the next error message
            messageStreamPtr_->rewind();

            throw errorException;
        }
        else
        {
            Perr<< endl << *this << endl
                << "\nFOAM aborting\n" << endl;
            printStack(Perr);
            ::abort();
        }
    }
}


Foam::Ostream& Foam::operator<<(Ostream& os, const IOerror& ioErr)
{
    if (!os.bad())
    {
        os  << endl
            << ioErr.title().c_str() << endl
            << ioErr.message().c_str() << endl << endl;

        os  << "file: " << ioErr.ioFileName().c_str();

        if (ioErr.ioStartLineNumber() >= 0 && ioErr.ioEndLineNumber() >= 0)
        {
            os  << " from line " << ioErr.ioStartLineNumber()
                << " to line " << ioErr.ioEndLineNumber() << '.';
        }
        else if (ioErr.ioStartLineNumber() >= 0)
        {
            os  << " at line " << ioErr.ioStartLineNumber() << '.';
        }

        if (IOerror::level >= 2 && ioErr.sourceFileLineNumber())
        {
            os  << endl << endl
                << "    From function " << ioErr.functionName().c_str() << endl
                << "    in file " << ioErr.sourceFileName().c_str()
                << " at line " << ioErr.sourceFileLineNumber() << '.';
        }
    }

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions

Foam::IOerror Foam::FatalIOError("--> FOAM FATAL IO ERROR: ");

// ************************************************************************* //
