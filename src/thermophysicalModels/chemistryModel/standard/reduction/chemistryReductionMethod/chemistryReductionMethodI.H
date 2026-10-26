/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2026 OpenFOAM Foundation
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

#include "chemistryReductionMethod.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::chemistryReductionMethod::nSpecie()
{
    return nSpecie_;
}


inline Foam::label
Foam::chemistryReductionMethod::nActiveSpecies() const
{
    return nActiveSpecies_;
}


inline const Foam::boolList&
Foam::chemistryReductionMethod::activeSpecies() const
{
    return activeSpecies_;
}


inline Foam::scalar
Foam::chemistryReductionMethod::tolerance() const
{
    return tolerance_;
}


inline bool Foam::chemistryReductionMethod::reactionDisabled
(
    const label i
) const
{
    return reactionsDisabled_[i];
}


// ************************************************************************* //
