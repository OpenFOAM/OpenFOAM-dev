/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "fvSpecificSource.H"
#include "polyCellSet.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvSpecificSource, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvSpecificSource::fvSpecificSource
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvSource(name, modelType, mesh, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvSpecificSource::~fvSpecificSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelUList Foam::fvSpecificSource::cells() const
{
    return polyCellSet(mesh()).cells();
}


Foam::label Foam::fvSpecificSource::nCells() const
{
    return mesh().nCells();
}


Foam::tmp<Foam::scalarField> Foam::fvSpecificSource::source
(
    const word& fieldName
) const
{
    tmp<volScalarField::Internal> tS(S(fieldName));

    return
        tS.isTmp()
      ? tmp<scalarField>(new scalarField(tS.ref(), true))
      : tmp<scalarField>(tS());
}


bool Foam::fvSpecificSource::movePoints()
{
    return true;
}


void Foam::fvSpecificSource::topoChange(const polyTopoChangeMap& map)
{}


void Foam::fvSpecificSource::mapMesh(const polyMeshMap& map)
{}


void Foam::fvSpecificSource::distribute(const polyDistributionMap& map)
{}


// ************************************************************************* //
