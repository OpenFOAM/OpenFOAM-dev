/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "perfectFluid.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::perfectFluid<Specie>::perfectFluid
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Specie(name, dict),
    R_(subDict.lookup<scalar>("R", dimSpecificHeatCapacity)),
    rho0_(subDict.lookup<scalar>("rho0", dimDensity))
{}


template<class Specie>
Foam::perfectFluid<Specie>::perfectFluid
(
    const word& name,
    const dictionary& dict
)
:
    perfectFluid(name, dict, dict.subDict("equationOfState"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::perfectFluid<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    writeEntry
    (
        os,
        "equationOfState",
        dictionary::entries("R", R_, "rho0", rho0_)
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const perfectFluid<Specie>& pf)
{
    pf.write(os);
    return os;
}


// ************************************************************************* //
