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

#include "icoPolynomial.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize>::icoPolynomial(Istream& is)
:
    Specie(is),
    rhoCoeffs_("rhoCoeffs<" + Foam::name(PolySize) + '>', is)
{
    rhoCoeffs_ *= this->W();
}


template<class Specie, int PolySize>
icoPolynomial<Specie, PolySize>::icoPolynomial(const dictionary& dict)
:
    Specie(dict),
    rhoCoeffs_
(
    dict.subDict("equationOfState").lookup
    (
        "rhoCoeffs<" + Foam::name(PolySize) + '>'
    )
)
{
    rhoCoeffs_ *= this->W();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
void icoPolynomial<Specie, PolySize>::write(Ostream& os) const
{
    Specie::write(os);

    dictionary dict("equationOfState");
    dict.add
    (
        word("rhoCoeffs<" + Foam::name(PolySize) + '>'),
        rhoCoeffs_/this->W()
    );

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Ostream& operator<<(Ostream& os, const icoPolynomial<Specie, PolySize>& ip)
{
    os  << static_cast<const Specie&>(ip) << tab
        << "rhoCoeffs<" << Foam::name(PolySize) << '>' << tab
        << ip.rhoCoeffs_/ip.W();

    os.check
    (
        "Ostream& operator<<"
        "(Ostream& os, const icoPolynomial<Specie, PolySize>& ip)"
    );

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
