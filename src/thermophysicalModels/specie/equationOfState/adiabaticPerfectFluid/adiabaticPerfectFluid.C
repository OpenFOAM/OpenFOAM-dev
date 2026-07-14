/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2026 OpenFOAM Foundation
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

#include "adiabaticPerfectFluid.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::adiabaticPerfectFluid<Specie>::adiabaticPerfectFluid
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Specie(name, dict),
    p0_(subDict.lookup<scalar>("p0", dimPressure)),
    rho0_(subDict.lookup<scalar>("rho0", dimDensity)),
    gamma_(subDict.lookup<scalar>("gamma", dimless)),
    B_(subDict.lookup<scalar>("B", dimPressure))
{}


template<class Specie>
Foam::adiabaticPerfectFluid<Specie>::adiabaticPerfectFluid
(
    const word& name,
    const dictionary& dict
)
:
    adiabaticPerfectFluid(name, dict, dict.subDict("equationOfState"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::adiabaticPerfectFluid<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    writeEntry
    (
        os,
        "equationOfState",
        dictionary::entries("p0", p0_, "rho0", rho0_, "gamma", gamma_, "B", B_)
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const adiabaticPerfectFluid<Specie>& pf
)
{
    pf.write(os);
    return os;
}


// ************************************************************************* //
