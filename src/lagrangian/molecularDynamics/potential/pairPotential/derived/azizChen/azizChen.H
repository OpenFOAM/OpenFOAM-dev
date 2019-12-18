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
    Foam::pairPotentials::azizChen

Description

    Reference:
    \verbatim
        Aziz, R. A. & Chen, H.H. (1977).
        An Accurate Intermolecular Potential for Argon.
        Journal of Chemical Physics, Vol. 67, No. 12, p. 5179.
    \endverbatim

SourceFiles
    azizChen.C

\*---------------------------------------------------------------------------*/

#ifndef azizChen_H
#define azizChen_H

#include "pairPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace pairPotentials
{

/*---------------------------------------------------------------------------*\
                           Class azizChen Declaration
\*---------------------------------------------------------------------------*/

class azizChen
:
    public pairPotential
{
    // Private Data

        dictionary azizChenCoeffs_;

        scalar epsilon_;
        scalar rm_;
        scalar A_;
        scalar alpha_;
        scalar C6_;
        scalar C8_;
        scalar C10_;
        scalar D_;
        scalar gamma_;


public:

    //- Runtime type information
    TypeName("azizChen");


    // Constructors

        //- Construct from components
        azizChen
        (
            const word& name,
            const dictionary& pairPotentialProperties
        );


    //- Destructor
    ~azizChen()
    {}


    // Member Functions

        scalar unscaledEnergy(const scalar r) const;

        //- Read dictionary
        bool read(const dictionary& pairPotentialProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
