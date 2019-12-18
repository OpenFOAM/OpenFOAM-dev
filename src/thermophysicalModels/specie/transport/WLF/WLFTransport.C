/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
Foam::scalar Foam::WLFTransport<Thermo>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").lookup<scalar>(coeffName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::WLFTransport<Thermo>::WLFTransport(const dictionary& dict)
:
    Thermo(dict),
    mu0_(readCoeff("mu0", dict)),
    Tr_(readCoeff("Tr", dict)),
    C1_(readCoeff("C1", dict)),
    C2_(readCoeff("C2", dict)),
    rPr_(1.0/dict.subDict("transport").lookup<scalar>("Pr"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::WLFTransport<Thermo>::write(Ostream& os) const
{
    os  << this->specie::name() << endl
        << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("mu0", mu0_);
    dict.add("Tr", Tr_);
    dict.add("C1", C1_);
    dict.add("C2", C2_);
    dict.add("Pr", 1.0/rPr_);

    os  << indent << dict.dictName() << dict
        << decrIndent << token::END_BLOCK << nl;
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
