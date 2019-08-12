/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "interpolatedWallDamping.H"
#include "phasePair.H"
#include "surfaceInterpolate.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallDampingModels
{
    defineTypeNameAndDebug(interpolated, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::wallDampingModels::interpolated::zeroNearWallCells
(
    const tmp<volScalarField>& tlimiter
) const
{
    if (zeroInNearWallCells_)
    {
        volScalarField& limiter = tlimiter.ref();

        const fvMesh& mesh(limiter.mesh());

        forAll(mesh.boundary(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                const labelUList& faceCells =
                    mesh.boundary()[patchi].faceCells();

                forAll(faceCells,facei)
                {
                    limiter[faceCells[facei]] = 0;
                }
            }
        }
    }

    return tlimiter;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDampingModels::interpolated::interpolated
(
    const dictionary& dict,
    const phasePair& pair
)
:
    wallDampingModel(dict, pair),
    zeroInNearWallCells_
    (
        dict.lookupOrDefault<Switch>("zeroInNearWallCells", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDampingModels::interpolated::~interpolated()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::wallDampingModels::interpolated::damp
(
    const tmp<volScalarField>& F
) const
{
    return zeroNearWallCells(limiter())*F;
}


Foam::tmp<Foam::volVectorField>
Foam::wallDampingModels::interpolated::damp
(
    const tmp<volVectorField>& F
) const
{
    return zeroNearWallCells(limiter())*F;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::wallDampingModels::interpolated::damp
(
    const tmp<surfaceScalarField>& Ff
) const
{
    return fvc::interpolate(zeroNearWallCells(limiter()))*Ff;
}


// ************************************************************************* //
