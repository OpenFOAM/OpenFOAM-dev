/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "hPowerThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::hPowerThermo<EquationOfState>::hPowerThermo
(
    const word& name,
    const dictionary& dict
)
:
    EquationOfState(name, dict),
    c0_(dict.subDict("thermodynamics").lookup<scalar>("C0")),
    n0_(dict.subDict("thermodynamics").lookup<scalar>("n0")),
    Tref_(dict.subDict("thermodynamics").lookup<scalar>("Tref")),
    hf_
    (
        dict
       .subDict("thermodynamics")
       .lookupBackwardsCompatible<scalar>({"hf", "Hf"})
    )
{}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const hPowerThermo<EquationOfState>& et
)
{
    et.write(os);
    return os;
}


// ************************************************************************* //
