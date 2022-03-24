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

#include "inverseDistanceDiffusivity.H"
#include "addToRunTimeSelectionTable.H"
#include "patchDistWave.H"
#include "wallPoint.H"
#include "HashSet.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inverseDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseDistanceDiffusivity::inverseDistanceDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    motionDiffusivity(mesh),
    patchNames_(mdData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseDistanceDiffusivity::~inverseDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::inverseDistanceDiffusivity::y() const
{
    const labelHashSet patchSet(mesh().boundaryMesh().patchSet(patchNames_));

    if (patchSet.size())
    {
        tmp<scalarField> tY(new scalarField(mesh().nCells()));
        patchDistWave::wave<wallPoint>(mesh(), patchSet, tY.ref(), false);
        return tY;
    }
    else
    {
        return tmp<scalarField>(new scalarField(mesh().nCells(), 1));
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::inverseDistanceDiffusivity::operator()() const
{
    volScalarField y_
    (
        IOobject
        (
            "y",
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimless,
        zeroGradientFvPatchScalarField::typeName
    );
    y_.primitiveFieldRef() = y();
    y_.correctBoundaryConditions();

    return surfaceScalarField::New
    (
        "faceDiffusivity",
        1.0/fvc::interpolate(y_)
    );
}


// ************************************************************************* //
