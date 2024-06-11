/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::OUprocess

Description
    Random Ornstein-Uhlenbeck process

SourceFiles
    OUprocess.C

\*---------------------------------------------------------------------------*/

#ifndef OUprocess_H
#define OUprocess_H

#include "complexFields.H"
#include "scalar.H"
#include "standardNormal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Kmesh;

/*---------------------------------------------------------------------------*\
                          Class OUprocess Declaration
\*---------------------------------------------------------------------------*/

class OUprocess
{
    // Private Data

        const distributions::standardNormal stdNormal_;

        const Kmesh& Kmesh_;
        mutable complexVectorField OUfield_;

        // Ornstein-Uhlenbeck process coefficients
        scalar alpha_;
        scalar sigma_;
        scalar kUpper_;
        scalar kLower_;
        scalar scale_;


    // Private Member Functions

        //- Return the Weiner process field
        complexVector WeinerProcess(const scalar deltaT) const;


public:

    // Constructors

        //- Construct from wavenumber mesh, timestep and coefficients dict
        OUprocess
        (
            const Kmesh& kmesh,
            const scalar deltaT,
            const dictionary&
        );


    // Member Functions

        //- Return the current random Ornstein-Uhlenbeck process field
        const complexVectorField& newField(const scalar deltaT) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
