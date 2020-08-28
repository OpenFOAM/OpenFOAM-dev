/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    Foam::functionObjects::logFiles

Description
    functionObject base class for creating, maintaining and writing log
    files e.g. integrated of averaged field data vs time.

See also
    Foam::functionObject
    Foam::functionObjects::writeFile

SourceFiles
    logFiles.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_logFiles_H
#define functionObjects_logFiles_H

#include "writeFile.H"
#include "OFstream.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class logFiles Declaration
\*---------------------------------------------------------------------------*/

class logFiles
:
    public writeFile
{

protected:

    // Protected data

        //- File names
        wordList names_;

        //- File pointer
        PtrList<OFstream> filePtrs_;


    // Protected Member Functions

        //- Create the output file
        virtual void createFiles();

        //- Reset the list of names from a wordList
        virtual void resetNames(const wordList& names);

        //- Reset the list of names to a single name entry
        virtual void resetName(const word& name);

        //- File header information
        virtual void writeFileHeader(const label i = 0) = 0;


public:

    // Constructors

        //- Construct from objectRegistry and prefix
        logFiles
        (
            const objectRegistry& obr,
            const word& prefix
        );

        //- Disallow default bitwise copy construction
        logFiles(const logFiles&) = delete;


    //- Destructor
    virtual ~logFiles();


    // Member Functions

        //- Return const access to the names
        const wordList& names() const;

        //- Return access to the files
        PtrList<OFstream>& files();

        //- Return access to the file (if only 1)
        OFstream& file();

        //- Return file 'i'
        OFstream& file(const label i);

        //- Write function
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const logFiles&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
