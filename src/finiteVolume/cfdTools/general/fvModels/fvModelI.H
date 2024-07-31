/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "fvModel.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline const Foam::dimensionSet& Foam::fvModel::sourceDims
(
    const dimensionSet& ds
)
{
    return ds;
}


inline const Foam::word& Foam::fvModel::fieldName()
{
    return word::null;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::fvModel::name() const
{
    return name_;
}


inline const Foam::word& Foam::fvModel::keyword() const
{
    return name_;
}


inline const Foam::fvMesh& Foam::fvModel::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary& Foam::fvModel::coeffs
(
    const dictionary& dict
) const
{
    return dict.optionalSubDict(modelType_ + "Coeffs");
}


// ************************************************************************* //
