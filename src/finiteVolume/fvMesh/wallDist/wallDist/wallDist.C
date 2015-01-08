/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "wallDist.H"
#include "wallPolyPatch.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallDist, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::wallDist::patchTypes(const labelHashSet& patchIDs) const
{
    wordList yTypes
    (
        mesh().boundary().size(),
        zeroGradientFvPatchField<Type>::typeName
    );

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        yTypes[iter.key()] = fixedValueFvPatchField<Type>::typeName;
    }

    return yTypes;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDist::wallDist(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::UpdateableMeshObject, wallDist>(mesh),
    pdm_
    (
        patchDistMethod::New
        (
            static_cast<const fvSchemes&>(mesh).subDict("wallDist"),
            mesh,
            mesh.boundaryMesh().findPatchIDs<wallPolyPatch>()
        )
    ),
    y_
    (
        IOobject
        (
            "yWall",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("yWall", dimLength, SMALL),
        patchTypes<scalar>(pdm_->patchIDs())
    ),
    n_(NULL)
{
    // Temporarily always construct n
    // until the demand-driven interface is complete
    n_ = tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "nWall",
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensionedVector("nWall", dimless, vector::zero),
            patchTypes<vector>(pdm_->patchIDs())
        )
    );

    const labelHashSet& patchIDs = pdm_->patchIDs();
    const fvPatchList& patches = mesh.boundary();

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        label patchi = iter.key();
        n_().boundaryField()[patchi] == patches[patchi].nf();
    }

    movePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDist::~wallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::wallDist::movePoints()
{
    if (pdm_->movePoints())
    {
        return pdm_->correct(y_, n_());
    }
    else
    {
        return false;
    }
}


void Foam::wallDist::updateMesh(const mapPolyMesh& mpm)
{
    pdm_->updateMesh(mpm);
    pdm_->correct(y_, n_());
}


// ************************************************************************* //
