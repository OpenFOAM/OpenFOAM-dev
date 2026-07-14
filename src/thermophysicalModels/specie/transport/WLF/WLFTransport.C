/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2026 OpenFOAM Foundation
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

#include "WLFTransport.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::WLFTransport<Thermo>::WLFTransport
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    Thermo(name, dict),
    mu0_(subDict.lookup<scalar>("mu0", dimDynamicViscosity)),
    Tr_(subDict.lookup<scalar>("Tr", dimTemperature)),
    C1_(subDict.lookup<scalar>("C1", dimless)),
    C2_(subDict.lookup<scalar>("C2", dimTemperature)),
    rPr_(1.0/subDict.lookup<scalar>("Pr", dimless))
{}


template<class Thermo>
Foam::WLFTransport<Thermo>::WLFTransport
(
    const word& name,
    const dictionary& dict
)
:
    WLFTransport(name, dict, dict.subDict("transport"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::WLFTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);

    writeEntry
    (
        os,
        "transport",
        dictionary::entries
        (
            "mu0", mu0_,
            "Tr", Tr_,
            "C1", C1_,
            "C2", C2_,
            "Pr", 1.0/rPr_
        )
    );
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const WLFTransport<Thermo>& wlft
)
{
    wlft.write(os);
    return os;
}


// ************************************************************************* //
