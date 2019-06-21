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

Class
    Foam::surfaceRegistry

Description
    Wraps the normal objectRegistry with a local instance for %surfaces.

SourceFiles
    surfaceRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceRegistry_H
#define surfaceRegistry_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class surfaceRegistry Declaration
\*---------------------------------------------------------------------------*/

class surfaceRegistry
:
    public objectRegistry
{
public:

    //- Runtime type information
    TypeName("surfaceRegistry");

        //- The prefix to local: %surfaces
        static const word prefix;

        //- The default surface name: %default
        static word defaultName;


    // Constructors

        //- Construct for the given objectRegistry and named surface
        surfaceRegistry(const objectRegistry&, const word& surfName = "");

        //- Disallow default bitwise copy construction
        surfaceRegistry(const surfaceRegistry&) = delete;


    //- Destructor
    virtual ~surfaceRegistry();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceRegistry&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
