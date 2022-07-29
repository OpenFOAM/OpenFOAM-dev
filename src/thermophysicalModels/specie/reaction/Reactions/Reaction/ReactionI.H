/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MulticomponentThermo>
inline Foam::scalar Foam::Reaction<MulticomponentThermo>::Tlow() const
{
    return Tlow_;
}


template<class MulticomponentThermo>
inline Foam::scalar Foam::Reaction<MulticomponentThermo>::Thigh() const
{
    return Thigh_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class MulticomponentThermo>
inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Reaction<MulticomponentThermo>& r
)
{
    r.write(os);
    return os;
}


// ************************************************************************* //
