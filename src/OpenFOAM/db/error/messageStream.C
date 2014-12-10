/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
#include "dictionary.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::messageStream::level(Foam::debug::debugSwitch("level", 2));


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::messageStream::messageStream
(
    const string& title,
    errorSeverity sev,
    const int maxErrors
)
:
    title_(title),
    severity_(sev),
    maxErrors_(maxErrors),
    errorCount_(0)
{}


Foam::messageStream::messageStream(const dictionary& dict)
:
    title_(dict.lookup("title")),
    severity_(FATAL),
    maxErrors_(0),
    errorCount_(0)
{}


Foam::OSstream& Foam::messageStream::masterStream(const label communicator)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** messageStream with comm:" << communicator
            << endl;
        error::printStack(Pout);
    }

    if (communicator == UPstream::worldComm)
    {
        return operator()();
    }
    else
    {
        return operator()(UPstream::master(communicator));
    }
}


Foam::OSstream& Foam::messageStream::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    OSstream& os = operator OSstream&();

    os  << endl
        << "    From function " << functionName << endl
        << "    in file " << sourceFileName
        << " at line " << sourceFileLineNumber << endl
        << "    ";

    return os;
}


Foam::OSstream& Foam::messageStream::operator()
(
    const string& functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber
)
{
    return operator()
    (
        functionName.c_str(),
        sourceFileName,
        sourceFileLineNumber
    );
}


Foam::OSstream& Foam::messageStream::operator()
(
    const char* functionName,
    const char* sourceFileName,
    const int sourceFileLineNumber,
    const string& ioFileName,
    const label ioStartLineNumber,
    const label ioEndLineNumber
)
{
    OSstream& os = operator OSstream&();

    os  << endl
        << "    From function " << functionName << endl
        << "    in file " << sourceFileName
        << " at line " << sourceFileLineNumber << endl
        << "    Reading " << ioFileName;

    if (ioStartLineNumber >= 0 && ioEndLineNumber >= 0)
    {
        os  << " from line " << ioStartLineNumber
            << " to line " << ioEndLineNumber;
    }
    else if (ioStartLineNumber >= 0)
    {
        os  << " at line " << ioStartLineNumber;
    }

    os << endl  << "    ";

    return os;
}


Foam::OSstream& Foam::messageStream::operator()
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


Foam::OSstream& Foam::messageStream::operator()
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


Foam::OSstream& Foam::messageStream::operator()(const bool output)
{
    if (output)
    {
        return operator()();
    }
    else
    {
        return Snull;
    }
}


Foam::messageStream::operator Foam::OSstream&()
{
    if (level)
    {
        bool collect = (severity_ == INFO || severity_ == WARNING);

        // Report the error
        if (!Pstream::master() && collect)
        {
            return Snull;
        }
        else
        {
            if (title().size())
            {
                if (Pstream::parRun() && !collect)
                {
                    Pout<< title().c_str();
                }
                else
                {
                    Sout<< title().c_str();
                }
            }

            if (maxErrors_)
            {
                errorCount_++;

                if (errorCount_ >= maxErrors_)
                {
                    FatalErrorIn("messageStream::operator OSstream&()")
                        << "Too many errors"
                        << abort(FatalError);
                }
            }

            if (Pstream::parRun() && !collect)
            {
                return Pout;
            }
            else
            {
                return Sout;
            }
        }
    }

    return Snull;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global messageStream definitions

Foam::messageStream Foam::SeriousError
(
    "--> FOAM Serious Error : ",
    messageStream::SERIOUS,
    100
);

Foam::messageStream Foam::Warning
(
    "--> FOAM Warning : ",
    messageStream::WARNING
);

Foam::messageStream Foam::Info("", messageStream::INFO);


// ************************************************************************* //
