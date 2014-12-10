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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionRate, class FallOffFunction>
inline Foam::FallOffReactionRate<ReactionRate, FallOffFunction>::
FallOffReactionRate
(
    const ReactionRate& k0,
    const ReactionRate& kInf,
    const FallOffFunction& F,
    const thirdBodyEfficiencies& tbes
)
:
    k0_(k0),
    kInf_(kInf),
    F_(F),
    thirdBodyEfficiencies_(tbes)
{}


template<class ReactionRate, class FallOffFunction>
inline Foam::FallOffReactionRate<ReactionRate, FallOffFunction>::
FallOffReactionRate
(
    const speciesTable& species,
    Istream& is
)
:
    k0_(species, is.readBegin("FallOffReactionRate(Istream&)")),
    kInf_(species, is),
    F_(is),
    thirdBodyEfficiencies_(species, is)
{
    is.readEnd("FallOffReactionRate(Istream&)");
}


template<class ReactionRate, class FallOffFunction>
inline Foam::FallOffReactionRate<ReactionRate, FallOffFunction>::
FallOffReactionRate
(
    const speciesTable& species,
    const dictionary& dict
)
:
    k0_(species, dict.subDict("k0")),
    kInf_(species, dict.subDict("kInf")),
    F_(dict.subDict("F")),
    thirdBodyEfficiencies_(species, dict.subDict("thirdBodyEfficiencies"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionRate, class FallOffFunction>
inline Foam::scalar
Foam::FallOffReactionRate<ReactionRate, FallOffFunction>::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c
) const
{
    scalar k0 = k0_(p, T, c);
    scalar kInf = kInf_(p, T, c);

    scalar Pr = k0*thirdBodyEfficiencies_.M(c)/kInf;

    return kInf*(Pr/(1 + Pr))*F_(T, Pr);
}


template<class ReactionRate, class FallOffFunction>
inline void Foam::FallOffReactionRate<ReactionRate, FallOffFunction>::write
(
    Ostream& os
) const
{
    os  << indent << "k0" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    k0_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "kInf" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    kInf_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "F" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    F_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;

    os  << indent << "thirdBodyEfficiencies" << nl;
    os  << indent << token::BEGIN_BLOCK << nl;
    os  << incrIndent;
    thirdBodyEfficiencies_.write(os);
    os  << decrIndent;
    os  << indent << token::END_BLOCK << nl;
}


template<class ReactionRate, class FallOffFunction>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const FallOffReactionRate<ReactionRate, FallOffFunction>& forr
)
{
    os  << token::BEGIN_LIST
        << forr.k0_ << token::SPACE
        << forr.kInf_ << token::SPACE
        << forr.F_ << token::SPACE
        << forr.thirdBodyEfficiencies_
        << token::END_LIST;
    return os;
}


// ************************************************************************* //
