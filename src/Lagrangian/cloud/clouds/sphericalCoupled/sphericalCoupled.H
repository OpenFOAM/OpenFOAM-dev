/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::clouds::sphericalCoupled

Description
    Base class for clouds of spherical particles which are coupled to a fluid

SourceFiles
    sphericalCoupled.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalCoupled_H
#define sphericalCoupled_H

#include "spherical.H"
#include "coupled.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class sphericalCoupled Declaration
\*---------------------------------------------------------------------------*/

class sphericalCoupled
{
public:

    // Public Static Data

        //- Run-time type information
        TypeName("sphericalCoupled");


    // Public Data

        //- Reynolds Number
        const CloudDerivedField<scalar>& Re;


    // Constructors

        //- Construct from a reference to the cloud and its base classes
        sphericalCoupled(const cloud&, const spherical&, const coupled&);


    //- Destructor
    virtual ~sphericalCoupled();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
