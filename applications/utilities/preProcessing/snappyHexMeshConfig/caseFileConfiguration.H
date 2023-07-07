/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::caseFileConfiguration

Description
    Base class for writing case files.

SourceFiles
    caseFileConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef caseFileConfiguration_H
#define caseFileConfiguration_H

#include "IOobject.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class caseFileConfiguration Declaration
\*---------------------------------------------------------------------------*/

class caseFileConfiguration
{

protected:

    // Protected Data

        //- IOobject to write file headers
        IOobject dict_;

        //- Output file stream to write file content
        OFstream os_;


    // Protected Functions

        //- Output stream to begin a new dictionary without name
        void beginDict(OFstream& os);

        //- Output stream to begin a new dictionary with name
        void beginDict(OFstream& os, const word& name);

        //- Output stream to end a dictionary
        void endDict(OFstream& os, bool newline=true);

        //- Output stream to begin a new list without name
        void beginList(OFstream& os);

        //- Output stream to begin a new list with name
        void beginList(OFstream& os, const word& name);

        //- Output stream to end a list
        void endList(OFstream& os, bool newline=true);


public:

    // Constructors

        //- Construct from components
        caseFileConfiguration
        (
            const fileName& name,
            const fileName& dir,
            const Time& time
        );

        //- Disallow default bitwise copy construction
        caseFileConfiguration(const caseFileConfiguration&) = delete;


    //- Destructor
    ~caseFileConfiguration();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const caseFileConfiguration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
