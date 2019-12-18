/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::adiabaticPerfectFluid<Specie>::adiabaticPerfectFluid
(
    const dictionary& dict
)
:
    Specie(dict),
    p0_(dict.subDict("equationOfState").lookup<scalar>("p0")),
    rho0_(dict.subDict("equationOfState").lookup<scalar>("rho0")),
    gamma_(dict.subDict("equationOfState").lookup<scalar>("gamma")),
    B_(dict.subDict("equationOfState").lookup<scalar>("B"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::adiabaticPerfectFluid<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("p0", p0_);
    dict.add("rho0", rho0_);
    dict.add("gamma", gamma_);
    dict.add("B", B_);

    os  << indent << dict.dictName() << dict;
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
