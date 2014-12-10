/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "constIsoSolidTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::constIsoSolidTransport<Thermo>::constIsoSolidTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    kappa_(readScalar(dict.subDict("transport").lookup("kappa")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::constIsoSolidTransport<Thermo>::constIsoSolidTransport::write
(
    Ostream& os
) const
{
    Thermo::write(os);

    dictionary dict("transport");
    dict.add("kappa", kappa_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const constIsoSolidTransport<Thermo>& ct
)
{
    operator<<(os, static_cast<const Thermo&>(ct));
    os << tab << ct.kappa_;

    os.check
    (
        "Ostream& operator<<(Ostream& os,const constIsoSolidTransport& ct)"
    );

    return os;
}


// ************************************************************************* //
