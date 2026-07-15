/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

#include "rhoTabulated.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::rhoTabulated<Specie>::rhoTabulated
(
    const word& name,
    const dictionary& dict
)
:
    Specie(name, dict),
    rho_
    (
        "rho",
        {dimensions::pressure, dimensions::temperature, dimensions::density},
        dict.subDict("equationOfState").subDict("rho")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::rhoTabulated<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    const delimitDictionary dlmt(os, "equationOfState"), dlmtRho(os, "rho");
    rho_.write
    (
        os,
        {dimensions::pressure, dimensions::temperature, dimensions::density}
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const rhoTabulated<Specie>& ip)
{
    ip.write(os);
    return os;
}


// ************************************************************************* //
