/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::SRIFallOffFunction::SRIFallOffFunction
(
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar d,
    const scalar e
)
:
    a_(a),
    b_(b),
    c_(c),
    d_(d),
    e_(e)
{}


inline Foam::SRIFallOffFunction::SRIFallOffFunction(Istream& is)
:
    a_(readScalar(is.readBegin("SRIFallOffFunction(Istream&)"))),
    b_(readScalar(is)),
    c_(readScalar(is)),
    d_(readScalar(is)),
    e_(readScalar(is))
{
    is.readEnd("SRIFallOffFunction(Istream&)");
}


inline Foam::SRIFallOffFunction::SRIFallOffFunction(const dictionary& dict)
:
    a_(readScalar(dict.lookup("a"))),
    b_(readScalar(dict.lookup("b"))),
    c_(readScalar(dict.lookup("c"))),
    d_(readScalar(dict.lookup("d"))),
    e_(readScalar(dict.lookup("e")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::SRIFallOffFunction::operator()
(
    const scalar T,
    const scalar Pr
) const
{
    scalar X = 1.0/(1.0 + sqr(log10(max(Pr, SMALL))));
    return d_*pow(a_*exp(-b_/T) + exp(-T/c_), X)*pow(T, e_);
}


inline void Foam::SRIFallOffFunction::write(Ostream& os) const
{
    os.writeKeyword("a") << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("b") << b_ << token::END_STATEMENT << nl;
    os.writeKeyword("c") << c_ << token::END_STATEMENT << nl;
    os.writeKeyword("d") << d_ << token::END_STATEMENT << nl;
    os.writeKeyword("e") << e_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::SRIFallOffFunction& srifof
)
{
    os  << token::BEGIN_LIST
        << srifof.a_
        << token::SPACE << srifof.b_
        << token::SPACE << srifof.c_
        << token::SPACE << srifof.d_
        << token::SPACE << srifof.e_
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
