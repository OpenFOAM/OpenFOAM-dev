/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

inline Foam::label
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
findFieldIndex
(
    const word& fieldName
) const
{
    return fieldIndices_[fieldName];
}


inline const Foam::List<Foam::dictionary>&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
output() const
{
    return output_;
}


inline const Foam::List<Foam::dictionary>&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
entries() const
{
    return entries_;
}


inline const Foam::List<Foam::scalar>&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
min() const
{
    return min_;
}


inline const Foam::List<Foam::label>&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
dim() const
{
    return dim_;
}


inline const Foam::List<Foam::scalar>&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
delta() const
{
    return delta_;
}


inline const Foam::List<Foam::scalar>&
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
max() const
{
    return max_;
}


inline Foam::word
Foam::radiationModels::absorptionEmissionModels::interpolationLookUpTable::
tableName() const
{
    return fileName_.name();
}


// ************************************************************************* //
