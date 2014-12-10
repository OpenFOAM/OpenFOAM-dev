/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Fundamental dimensioned constants

\*---------------------------------------------------------------------------*/

#ifndef fundamentalConstants_H
#define fundamentalConstants_H

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace universal
{
    //- Speed of light in a vacuum
    extern const dimensionedScalar c;

    //- Newtonian constant of gravitation
    extern const dimensionedScalar G;

    //- Planck constant
    extern const dimensionedScalar h;
}

namespace electromagnetic
{
    //- Elementary charge
    extern const dimensionedScalar e;
}

namespace atomic
{
    //- Electron mass
    extern const dimensionedScalar me;

    //- Proton mass
    extern const dimensionedScalar mp;
}

namespace physicoChemical
{
    //- Atomic mass unit
    extern const dimensionedScalar mu;

    //- Avagadro number
    extern const dimensionedScalar NA;

    //- Boltzmann constant
    extern const dimensionedScalar k;
}

namespace standard
{
    //- Standard pressure
    extern const dimensionedScalar Pstd;

    //- Standard temperature
    extern const dimensionedScalar Tstd;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace constant
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
