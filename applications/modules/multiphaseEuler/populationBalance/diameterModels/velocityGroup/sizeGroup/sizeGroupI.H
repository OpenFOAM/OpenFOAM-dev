/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

#include "sizeGroup.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::diameterModels::sizeGroup::i() const
{
    return i_;
}


inline const Foam::phaseModel&
Foam::diameterModels::sizeGroup::phase() const
{
    return group().phase();
}


inline const Foam::diameterModels::velocityGroup&
Foam::diameterModels::sizeGroup::group() const
{
    return group_;
}


inline const Foam::dimensionedScalar&
Foam::diameterModels::sizeGroup::dSph() const
{
    return dSph_;
}


inline const Foam::dimensionedScalar&
Foam::diameterModels::sizeGroup::x() const
{
    return x_;
}


inline Foam::autoPtr<Foam::diameterModels::shapeModel>&
Foam::diameterModels::sizeGroup::shapeModelPtr()
{
    return shapeModel_;
}


inline const Foam::autoPtr<Foam::diameterModels::shapeModel>&
Foam::diameterModels::sizeGroup::shapeModelPtr() const
{
    return shapeModel_;
}

// ************************************************************************* //
