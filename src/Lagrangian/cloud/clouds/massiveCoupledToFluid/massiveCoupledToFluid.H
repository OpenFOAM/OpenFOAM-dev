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
    Foam::clouds::massiveCoupledToFluid

Description
    Base class for clouds of massive particles which are coupled to a fluid

SourceFiles
    massiveCoupledToFluid.C

\*---------------------------------------------------------------------------*/

#ifndef massiveCoupledToFluid_H
#define massiveCoupledToFluid_H

#include "massive.H"
#include "coupledToFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class massiveCoupledToFluid Declaration
\*---------------------------------------------------------------------------*/

class massiveCoupledToFluid
{
public:

    // Public Static Data

        //- Run-time type information
        TypeName("massiveCoupledToFluid");


    // Public Data

        //- Cloud/carrier density ratio
        const CloudDerivedField<scalar>& rhoByRhoc;


    // Constructors

        //- Construct from a reference to the cloud and its base classes
        massiveCoupledToFluid
        (
            const cloud&,
            const massive&,
            const coupledToFluid&
        );


    //- Destructor
    virtual ~massiveCoupledToFluid();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
