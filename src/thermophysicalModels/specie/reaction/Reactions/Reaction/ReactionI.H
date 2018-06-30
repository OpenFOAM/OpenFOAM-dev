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

#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
inline const word& Reaction<ReactionThermo>::name() const
{
    return name_;
}


template<class ReactionThermo>
inline scalar Reaction<ReactionThermo>::Tlow() const
{
    return Tlow_;
}


template<class ReactionThermo>
inline scalar Reaction<ReactionThermo>::Thigh() const
{
    return Thigh_;
}


template<class ReactionThermo>
inline const List<typename Reaction<ReactionThermo>::specieCoeffs>&
Reaction<ReactionThermo>::lhs() const
{
    return lhs_;
}


template<class ReactionThermo>
inline const List<typename Reaction<ReactionThermo>::specieCoeffs>&
Reaction<ReactionThermo>::rhs() const
{
    return rhs_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class ReactionThermo>
inline Ostream& operator<<(Ostream& os, const Reaction<ReactionThermo>& r)
{
    r.write(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
