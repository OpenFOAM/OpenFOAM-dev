/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "movingMappedWallVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "mappedFvPatchBaseBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::movingMappedWallVelocityFvPatchVectorField::
movingMappedWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict)
{}


Foam::movingMappedWallVelocityFvPatchVectorField::
movingMappedWallVelocityFvPatchVectorField
(
    const movingMappedWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper, false)
{
    mapper(*this, ptf, vector::zero);
}


Foam::movingMappedWallVelocityFvPatchVectorField::
movingMappedWallVelocityFvPatchVectorField
(
    const movingMappedWallVelocityFvPatchVectorField& mwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(mwvpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::movingMappedWallVelocityFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    mapper(*this, ptf, vector::zero);
}


void Foam::movingMappedWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatch& fvp = patch();

    const mappedFvPatchBaseBase& mapper =
        mappedFvPatchBaseBase::getMap(patch());
    const fvMesh& nbrMesh = mapper.nbrMesh();
    const fvPatch& nbrFvp = mapper.nbrFvPatch();

    if (nbrMesh.moving())
    {
        // Calculate the normal velocity on this mesh from the mesh flux
        const volVectorField& U =
            static_cast<const volVectorField&>(internalField());
        tmp<scalarField> Un =
            U.mesh().moving()
          ? fvc::meshPhi(U, fvp.index())/(fvp.magSf() + vSmall)
          : tmp<scalarField>(new scalarField(fvp.size(), Zero));

        // Calculate and map the mesh velocity from the neighbour
        vectorField nbrCf(nbrFvp.size()), nbrCf0(nbrFvp.size());
        forAll(nbrCf0, nbrPatchFacei)
        {
            const label nbrPolyFacei =
                nbrMesh.polyFacesBf()[nbrFvp.index()][nbrPatchFacei];
            nbrCf[nbrPatchFacei] =
                nbrMesh.faceCentres()[nbrPolyFacei];
            nbrCf0[nbrPatchFacei] =
                nbrMesh.faces()[nbrPolyFacei].centre(nbrMesh.oldPoints());
        }
        const vectorField nbrUp
        (
            (nbrCf - nbrCf0)/db().time().deltaTValue()
        );
        const vectorField Up(mapper.fromNeighbour(nbrUp));

        // Set the velocity equal to the normal component from this mesh plus
        // the tangential component from the neighbouring mesh. That way the
        // flux is exactly preserved but the neighbour can impart shear. The
        // normal components should be the same anyway, otherwise the mapped
        // patches are moving apart. Using the normal component from this mesh
        // just prevents mapping inaccuracies from affecting the flux.
        const vectorField n(fvp.nf());
        vectorField::operator=(Up + n*(Un - (n & Up)));
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::movingMappedWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        movingMappedWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
