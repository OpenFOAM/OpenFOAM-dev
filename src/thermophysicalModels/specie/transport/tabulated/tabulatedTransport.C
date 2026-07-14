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

#include "tabulatedTransport.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Thermo(name, dict),
    mu_
    (
        "mu",
        {dimPressure, dimTemperature, dimDynamicViscosity},
        subDict.subDict("mu")
    ),
    kappa_
    (
        "kappa",
        {dimPressure, dimTemperature, dimThermalConductivity},
        subDict.subDict("kappa")
    )
{}


template<class Thermo>
Foam::tabulatedTransport<Thermo>::tabulatedTransport
(
    const word& name,
    const dictionary& dict
)
:
    tabulatedTransport(name, dict, dict.subDict("transport"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::tabulatedTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);

    const delimitDictionary delimit(os, "transport");
    {
        const delimitDictionary delimitMu(os, "mu");
        mu_.write(os, {dimPressure, dimTemperature, dimDynamicViscosity});
    }
    {
        const delimitDictionary delimitKappa(os, "kappa");
        kappa_.write(os, {dimPressure, dimTemperature, dimThermalConductivity});
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tabulatedTransport<Thermo>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
