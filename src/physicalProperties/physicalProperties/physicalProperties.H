/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::physicalProperties

Description
    A base class for physical properties.

SourceFiles
    physicalProperties.C

\*---------------------------------------------------------------------------*/

#ifndef physicalProperties_H
#define physicalProperties_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class physicalProperties Declaration
\*---------------------------------------------------------------------------*/

class physicalProperties
:
    public IOdictionary
{

public:

    //- Runtime type information
    TypeName("physicalProperties");


    // Constructors

        //- Construct from components
        physicalProperties(const objectRegistry& obr, const word& group);

        //- Disallow default bitwise copy construction
        physicalProperties(const physicalProperties&) = delete;


    //- Destructor
    virtual ~physicalProperties()
    {}


    // Member Functions

        //- Find the physicalProperties dictionary and return the IOobject
        static IOobject findModelDict
        (
            const objectRegistry& obr,
            const word& group,
            bool registerObject = false
        );

        //- Read physicalProperties dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const physicalProperties&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
