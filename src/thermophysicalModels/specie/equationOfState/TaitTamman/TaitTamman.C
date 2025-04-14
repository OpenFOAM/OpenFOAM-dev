/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "TaitTamman.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::TaitTamman<Specie>::TaitTamman
(
    const word& name,
    const dictionary& dict
)
:
    Specie(name, dict),
    T0_(dict.subDict("equationOfState").lookup<scalar>("T0")),
    b1_(dict.subDict("equationOfState").lookup<scalar>("b1")),
    b2_(dict.subDict("equationOfState").lookup<scalar>("b2")),
    b3_(dict.subDict("equationOfState").lookup<scalar>("b3")),
    b4_(dict.subDict("equationOfState").lookup<scalar>("b4"))
{}

template<class Specie>
Foam::scalar Foam::TaitTamman<Specie>::C_ = 0.0894;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::TaitTamman<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("T0", T0_);
    dict.add("b1", b1_);
    dict.add("b2", b2_);
    dict.add("b3", b3_);
    dict.add("b4", b4_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TaitTamman<Specie>& pf
)
{
    pf.write(os);
    return os;
}


// ************************************************************************* //
