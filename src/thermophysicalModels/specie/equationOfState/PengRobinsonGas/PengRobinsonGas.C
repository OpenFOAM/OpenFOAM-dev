/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "PengRobinsonGas.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::PengRobinsonGas<Specie>::PengRobinsonGas(Istream& is)
:
    Specie(is),
    Tc_(readScalar(is)),
    Vc_(readScalar(is)),
    Zc_(readScalar(is)),
    Pc_(readScalar(is)),
    omega_(readScalar(is))
{
    is.check("PengRobinsonGas<Specie>::PengRobinsonGas(Istream& is)");
}


template<class Specie>
Foam::PengRobinsonGas<Specie>::PengRobinsonGas
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(readScalar(dict.subDict("equationOfState").lookup("Tc"))),
    Vc_(readScalar(dict.subDict("equationOfState").lookup("Vc"))),
    Zc_(1.0),
    Pc_(readScalar(dict.subDict("equationOfState").lookup("Pc"))),
    omega_(readScalar(dict.subDict("equationOfState").lookup("omega")))
{
    Zc_ = Pc_*Vc_/(specie::RR*Tc_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::PengRobinsonGas<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PengRobinsonGas<Specie>& pg
)
{
    os  << static_cast<const Specie&>(pg)
        << token::SPACE << pg.Tc_
        << token::SPACE << pg.Vc_
        << token::SPACE << pg.Zc_
        << token::SPACE << pg.Pc_
        << token::SPACE << pg.omega_;

    os.check
    (
        "Ostream& operator<<(Ostream& os, const PengRobinsonGas<Specie>& st)"
    );
    return os;
}


// ************************************************************************* //
