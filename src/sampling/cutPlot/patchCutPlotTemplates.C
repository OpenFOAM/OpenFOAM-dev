/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

#include "patchCutPlot.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class PatchType>
Foam::List<Foam::patchCutPlot::weight>
Foam::patchCutPlot::calcWeights
(
    const PatchType& p,
    const scalarField& localPointXs,
    const scalarField& cutXs,
    const bool interpolate,
    const bool normalise
)
{
    return
        calcWeights
        (
            p.localFaces(),
            p.faceAreas(),
            p.faceNormals(),
            p.localPoints(),
            localPointXs,
            cutXs,
            interpolate,
            normalise
        );
}


template<class PatchType>
Foam::tmp<Foam::scalarField> Foam::patchCutPlot::calcCutXs
(
    const PatchType& p,
    const scalarField& localPointXs,
    const bool interpolate,
    const label nCuts,
    const label nIter,
    const bool debug,
    const word& functionName,
    const polyMesh& functionMesh,
    const setWriter& functionFormatter
)
{
    return
        calcCutXs
        (
            p.localFaces(),
            p.faceAreas(),
            p.faceNormals(),
            p.localPoints(),
            localPointXs,
            interpolate,
            nCuts,
            nIter,
            debug,
            functionName,
            functionMesh,
            functionFormatter
        );
}


// ************************************************************************* //
