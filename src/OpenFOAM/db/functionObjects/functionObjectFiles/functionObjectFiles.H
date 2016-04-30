/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::functionObjectFiles

Description
    Base class for output file data handling

See Also
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    functionObjectFiles.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectFiles_H
#define functionObjectFiles_H

#include "functionObjectFile.H"
#include "OFstream.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class functionObjectFiles Declaration
\*---------------------------------------------------------------------------*/

class functionObjectFiles
:
    public functionObjectFile
{
    // Private data

        //- File names
        wordList names_;

        //- File pointer
        PtrList<OFstream> filePtrs_;


protected:

    // Protected Member Functions

        //- Create the output file
        virtual void createFiles();

        //- Write function
        virtual void write();

        //- Reset the list of names from a wordList
        virtual void resetNames(const wordList& names);

        //- Reset the list of names to a single name entry
        virtual void resetName(const word& name);

        //- Disallow default bitwise copy construct
        functionObjectFiles(const functionObjectFiles&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectFiles&);


public:

    // Constructors

        //- Construct from objectRegistry
        functionObjectFiles(const objectRegistry& obr, const word& prefix);

        //- Construct from components
        functionObjectFiles
        (
            const objectRegistry& obr,
            const word& prefix,
            const word& name
        );

        //- Construct from components
        functionObjectFiles
        (
            const objectRegistry& obr,
            const word& prefix,
            const wordList& names
        );


    //- Destructor
    virtual ~functionObjectFiles();


    // Member Functions

        //- Return const access to the names
        const wordList& names() const;

        //- Return access to the file (if only 1)
        OFstream& file();

        //- Return access to the files
        PtrList<OFstream>& files();

        //- Return file 'i'
        OFstream& file(const label i);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
