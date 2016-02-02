/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "logPolynomialTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::logPolynomialTransport<Thermo, PolySize>::logPolynomialTransport
(
    Istream& is
)
:
    Thermo(is),
    muCoeffs_("muLogCoeffs<" + Foam::name(PolySize) + '>', is),
    kappaCoeffs_("kappaLogCoeffs<" + Foam::name(PolySize) + '>', is)
{
    muCoeffs_ *= this->W();
    kappaCoeffs_ *= this->W();
}


template<class Thermo, int PolySize>
Foam::logPolynomialTransport<Thermo, PolySize>::logPolynomialTransport
(
    const dictionary& dict
)
:
    Thermo(dict),
    muCoeffs_
    (
        dict.subDict("transport").lookup
        (
            "muLogCoeffs<" + Foam::name(PolySize) + '>'
        )
    ),
    kappaCoeffs_
    (
        dict.subDict("transport").lookup
        (
            "kappaLogCoeffs<" + Foam::name(PolySize) + '>'
        )
    )
{
    muCoeffs_ *= this->W();
    kappaCoeffs_ *= this->W();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
void Foam::logPolynomialTransport<Thermo, PolySize>::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add
    (
        word("muLogCoeffs<" + Foam::name(PolySize) + '>'),
        muCoeffs_/this->W()
    );
    dict.add
    (
        word("kappaLogCoeffs<" + Foam::name(PolySize) + '>'),
        kappaCoeffs_/this->W()
    );
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const logPolynomialTransport<Thermo, PolySize>& pt
)
{
    os  << static_cast<const Thermo&>(pt) << tab
        << "muLogCoeffs<" << Foam::name(PolySize) << '>' << tab
        << pt.muCoeffs_/pt.W() << tab
        << "kappaLogCoeffs<" << Foam::name(PolySize) << '>' << tab
        << pt.kappaCoeffs_/pt.W();

    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const logPolynomialTransport<Thermo, PolySize>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
