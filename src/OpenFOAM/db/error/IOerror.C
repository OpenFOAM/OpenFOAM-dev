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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "OStringStream.H"
#include "fileName.H"
#include "dictionary.H"
#include "jobInfo.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::IOerrorLocation::IOerrorLocation()
:
    ioFileName_("unknown"),
    ioLineNumber_(-1),
    ioGlobal_(false)
{}


Foam::IOerrorLocation::IOerrorLocation
(
    const string& ioFileName,
    const label ioLineNumber,
    const bool ioGlobal
)
:
    ioFileName_(ioFileName),
    ioLineNumber_(ioLineNumber),
    ioGlobal_(ioGlobal)
{}


Foam::IOerrorLocation::IOerrorLocation(const IOstream& ios)
:
    ioFileName_(ios.name()),
    ioLineNumber_(ios.lineNumber()),
    ioGlobal_(ios.global())
{}


Foam::IOerrorLocation::IOerrorLocation(const dictionary& dict)
:
    ioFileName_(dict.name()),
    ioLineNumber_(dict.endLineNumber()),
    ioGlobal_(dict.global())
{}


Foam::IOerror::IOerror(const string& title)
:
    error(title),
    IOerrorLocation()
{}


Foam::OSstream& Foam::IOerror::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const IOerrorLocation& location
)
{
    error::operator()(functionName, sourceFileName, sourceFileLineNumber);

    IOerrorLocation::operator=(location);

    return operator OSstream&();
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
    errDict.add("ioLineNumber", ioLineNumber());

    return errDict;
}


void Foam::IOerror::exit(const int errNo)
{
    if (IOerror::level <= 0)
    {
        if (Pstream::parRun())
        {
            Pstream::exit(errNo);
        }
        else
        {
            ::exit(errNo);
        }
    }

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
        if (ioGlobal())
        {
            if (Pstream::master())
            {
                Serr<< endl << *this << endl
                    << "\nFOAM parallel run exiting\n" << endl;
            }
        }
        else
        {
            Perr<< endl << *this << endl
                << "\nFOAM parallel run exiting\n" << endl;
        }

        Pstream::exit(errNo);
    }
    else
    {
        if (throwExceptions_)
        {
            // Make a copy of the error to throw
            IOerror errorException(*this);

            // Rewind the message buffer for the next error message
            messageStream_.rewind();

            throw errorException;
        }
        else
        {
            Serr<< endl << *this << endl
                << "\nFOAM exiting\n" << endl;
            ::exit(errNo);
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
        if (ioGlobal())
        {
            if (Pstream::master())
            {
                Serr<< endl << *this << endl
                    << "\nFOAM parallel run aborting\n" << endl;
                printStack(Perr);
            }
        }
        else
        {
            Perr<< endl << *this << endl
                << "\nFOAM parallel run aborting\n" << endl;
            printStack(Perr);
        }

        Pstream::abort();
    }
    else
    {
        if (throwExceptions_)
        {
            // Make a copy of the error to throw
            IOerror errorException(*this);

            // Rewind the message buffer for the next error message
            messageStream_.rewind();

            throw errorException;
        }
        else
        {
            Serr<< endl << *this << endl
                << "\nFOAM aborting\n" << endl;
            printStack(Serr);
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

        if (ioErr.ioLineNumber() >= 0)
        {
            os  << " at line " << ioErr.ioLineNumber() << '.';
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
