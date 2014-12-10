/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "linear.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::linear<Specie>::linear(Istream& is)
:
    Specie(is),
    psi_(readScalar(is)),
    rho0_(readScalar(is))
{
    is.check("linear<Specie>::linear(Istream& is)");
}


template<class Specie>
Foam::linear<Specie>::linear(const dictionary& dict)
:
    Specie(dict),
    psi_(readScalar(dict.subDict("equationOfState").lookup("psi"))),
    rho0_(readScalar(dict.subDict("equationOfState").lookup("rho0")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::linear<Specie>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add("psi", psi_);
    dict.add("rho0", rho0_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<(Ostream& os, const linear<Specie>& pf)
{
    os  << static_cast<const Specie&>(pf)
        << token::SPACE << pf.psi_
        << token::SPACE << pf.rho0_;

    os.check("Ostream& operator<<(Ostream&, const linear<Specie>&)");
    return os;
}


// ************************************************************************* //
