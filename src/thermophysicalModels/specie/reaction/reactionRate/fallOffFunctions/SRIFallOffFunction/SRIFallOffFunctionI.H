/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    scalar X = 1.0/(1 + sqr(log10(max(Pr, small))));
    return d_*pow(a_*exp(-b_/T) + exp(-T/c_), X)*pow(T, e_);
}


inline Foam::scalar Foam::SRIFallOffFunction::ddT
(
    const scalar Pr,
    const scalar F,
    const scalar dPrdT,
    const scalar T
) const
{
    scalar X = 1.0/(1 + sqr(log10(max(Pr, small))));
    scalar dXdPr = -X*X*2*log10(Pr)/Pr/log(10.0);
    return
        F
       *(
            e_/T
          + X
           *(a_*b_*exp(-b_/T)/sqr(T) - exp(-T/c_)/c_)
           /(a_*exp(-b_/T) + exp(-T/c_))
          + dXdPr*dPrdT*log(a_*exp(-b_/T) + exp(-T/c_))
        );
}


inline Foam::scalar Foam::SRIFallOffFunction::ddc
(
    const scalar Pr,
    const scalar F,
    const scalar dPrdc,
    const scalar T
) const
{
    scalar X = 1.0/(1 + sqr(log10(max(Pr, small))));
    scalar dXdPr = -X*X*2*log10(Pr)/Pr/log(10.0);
    scalar dXdc = dXdPr*dPrdc;
    return F*dXdc*log(a_*exp(-b_/T) + exp(-T/c_));
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
    srifof.write(os);
    return os;
}


// ************************************************************************* //
