/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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

#include "wallDampingModel.H"
#include "surfaceInterpolate.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallDampingModel, 0);
    defineRunTimeSelectionTable(wallDampingModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDampingModel::wallDampingModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    wallDependentModel(interface.mesh()),
    interface_
    (
        interface.modelCast<wallDampingModel, dispersedPhaseInterface>()
    ),
    Cd_("Cd", dimless, dict),
    zeroWallDist_("zeroWallDist", dimLength, dict, 0),
    zeroInNearWallCells_
    (
        dict.lookupOrDefault<Switch>("zeroInNearWallCells", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDampingModel::~wallDampingModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::wallDampingModel::damping() const
{
    tmp<volScalarField> tlimiter(limiter().ptr());

    if (zeroInNearWallCells_)
    {
        volScalarField& limiter = tlimiter.ref();

        const fvBoundaryMesh& bMesh = limiter.mesh().boundary();

        forAll(bMesh, patchi)
        {
            if (isA<wallFvPatch>(bMesh[patchi]))
            {
                const labelUList& faceCells = bMesh[patchi].faceCells();

                forAll(faceCells, facei)
                {
                    limiter[faceCells[facei]] = 0;
                }
            }
        }

        return tlimiter;
    }
    else
    {
        return tlimiter;
    }
}


Foam::tmp<Foam::surfaceScalarField>
Foam::wallDampingModel::dampingf() const
{
    return fvc::interpolate(damping());
}


// ************************************************************************* //
