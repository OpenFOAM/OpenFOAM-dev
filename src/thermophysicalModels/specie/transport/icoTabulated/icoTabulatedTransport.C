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

#include "icoTabulatedTransport.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::icoTabulatedTransport<Thermo>::icoTabulatedTransport
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
        {dimensions::temperature, dimensions::dynamicViscosity},
        subDict.subDict("mu")
    ),
    kappa_
    (
        "kappa",
        {dimensions::temperature, dimensions::thermalConductivity},
        subDict.subDict("kappa")
    )
{}


template<class Thermo>
Foam::icoTabulatedTransport<Thermo>::icoTabulatedTransport
(
    const word& name,
    const dictionary& dict
)
:
    icoTabulatedTransport(name, dict, dict.subDict("transport"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::icoTabulatedTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);

    const delimitDictionary delimit(os, "transport");
    {
        const delimitDictionary delimitMu(os, "mu");
        mu_.write(os, {dimensions::temperature, dimensions::dynamicViscosity});
    }
    {
        const delimitDictionary delimitKappa(os, "kappa");
        kappa_.write
        (
            os,
            {dimensions::temperature, dimensions::thermalConductivity}
        );
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const icoTabulatedTransport<Thermo>& pt
)
{
    pt.write(os);
    return os;
}


// ************************************************************************* //
