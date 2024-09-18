/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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
    Foam::functionObjects::regionFunctionObject

Description
    Specialisation of Foam::functionObject for a region

See also
    Foam::functionObject
    Foam::functionObjects::objectRegistryFunctionObject

SourceFiles
    regionFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_regionFunctionObject_H
#define functionObjects_regionFunctionObject_H

#include "objectRegistryFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class regionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class regionFunctionObject
:
    public objectRegistryFunctionObject
{
public:

    //- Runtime type information
    TypeName("regionFunctionObject");


    // Constructors

        //- Construct from Time and dictionary.
        //  The region objectRegistry is looked-up runTime with the name
        //  looked-up from the dictionary (defaults to polyMesh::defaultRegion)
        regionFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from the region objectRegistry
        regionFunctionObject
        (
            const word& name,
            const objectRegistry& obr
        );

        //- Disallow default bitwise copy construction
        regionFunctionObject(const regionFunctionObject&) = delete;


    //- Destructor
    virtual ~regionFunctionObject();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const regionFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
