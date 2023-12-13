/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Namespace
    Foam::constant::physicoChemical

Description
    Physico-chemical constants

\*---------------------------------------------------------------------------*/

#ifndef physicoChemicalConstants_H
#define physicoChemicalConstants_H

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{
namespace physicoChemical
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Group name for physico-chemical constants
    extern const char* const group;

    //- Universal gas constant: default SI units: [J/kmol/K]
    extern const dimensionedScalar RR;

    //- Faraday constant: default SI units: [C/kmol]
    extern const dimensionedScalar F;

    //- Stefan-Boltzmann constant: default SI units: [W/m^2/K^4]
    extern const dimensionedScalar sigma;

    //- Wien displacement law constant: default SI units: [m K]
    extern const dimensionedScalar b;

    //- First radiation constant: default SI units: [W/m^2]
    extern const dimensionedScalar c1;

    //- Second radiation constant: default SI units: [m K]
    extern const dimensionedScalar c2;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace physicoChemical
} // End namespace constant
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
