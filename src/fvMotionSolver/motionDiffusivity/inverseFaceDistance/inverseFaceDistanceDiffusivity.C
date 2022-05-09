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

#include "inverseFaceDistanceDiffusivity.H"
#include "surfaceFields.H"
#include "fvPatchDistWave.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseFaceDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inverseFaceDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseFaceDistanceDiffusivity::inverseFaceDistanceDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    motionDiffusivity(mesh),
    patchNames_(mdData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseFaceDistanceDiffusivity::~inverseFaceDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::inverseFaceDistanceDiffusivity::operator()() const
{
    surfaceScalarField y
    (
        IOobject
        (
            "y",
            mesh().time().timeName(),
            mesh()
        ),
        mesh(),
        dimensionedScalar(dimless, 1)
    );

    const labelHashSet patchIDs(mesh().boundaryMesh().patchSet(patchNames_));

    if (patchNames_.size())
    {
        fvPatchDistWave::calculate
        (
            mesh(),
            patchIDs,
            -vGreat,
            y
        );

        // Use cell distance on faces that are part of the patch set. This
        // avoids divide-by-zero issues.
        forAllConstIter(labelHashSet, patchIDs, iter)
        {
            const label patchi = iter.key();

            const labelUList& patchCells =
                mesh().boundary()[patchi].faceCells();

            forAll(patchCells, patchFacei)
            {
                y.boundaryFieldRef()[patchi][patchFacei] =
                    mag
                    (
                       mesh().Cf().boundaryField()[patchi][patchFacei]
                     - mesh().C()[patchCells[patchFacei]]
                    );
            }
        }
    }

    return surfaceScalarField::New("faceDiffusivity", 1/y);
}


// ************************************************************************* //
