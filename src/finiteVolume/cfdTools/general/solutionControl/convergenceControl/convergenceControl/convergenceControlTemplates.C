/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ResidualData>
Foam::label Foam::convergenceControl::residualControlIndex
(
    const word& fieldName,
    const List<ResidualData>& residualControl,
    const bool useRegEx
)
{
    forAll(residualControl, i)
    {
        if (useRegEx && residualControl[i].name.match(fieldName))
        {
            return i;
        }
        else if (residualControl[i].name == fieldName)
        {
            return i;
        }
    }

    return -1;
}


template<class Type>
void Foam::convergenceControl::getInitialTypeResiduals
(
    const fvMesh& mesh,
    const word& fieldName,
    const label solvei,
    ITstream& data,
    scalar& r0,
    scalar& r
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (mesh.foundObject<fieldType>(fieldName))
    {
        const List<SolverPerformance<Type>> sp(data);
        r0 = cmptMax(sp[0].initialResidual());
        r = cmptMax(sp[solvei].initialResidual());
    }
}


// ************************************************************************* //
