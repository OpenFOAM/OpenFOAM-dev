/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::logFile

Description
    functionObject support for writing log files

    e.g. residuals or integral field values vs time.

See also
    Foam::functionObject
    Foam::functionObjects::writeFile

SourceFiles
    logFile.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_logFile_H
#define functionObjects_logFile_H

#include "writeFile.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class logFile Declaration
\*---------------------------------------------------------------------------*/

class logFile
:
    public writeFile
{
    // Private data

        //- File name
        word name_;

        //- File stream
        OFstream file_;


    // Private member functions

        fileName filePathName() const;


public:

    // Constructors

        //- Construct from objectRegistry, prefix and file name
        logFile
        (
            const objectRegistry& obr,
            const word& prefix,
            const word& name
        );

        //- Disallow default bitwise copy construction
        logFile(const logFile&) = delete;


    //- Destructor
    ~logFile();


    // Member Functions

        //- Return the log file name
        const word& name() const;

        //- Write the time and column headers
        void writeTimeColumnHeaders(const wordList& titles);

        //- Write a single column value
        template<class T>
        void write(const T& t)
        {
            file_ << tab << setw(charWidth()) << t << endl;
        }

        //- Write the column values
        template<class T, class... Args>
        void write(const T& t, const Args& ... args)
        {
            file_ << tab << setw(charWidth()) << t;
            write(args...);
        }

        //- Write the time and column values
        template<class... Args>
        void writeTime(const Args& ... args)
        {
            writeFile::writeTime(file_);
            write(args...);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const logFile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
