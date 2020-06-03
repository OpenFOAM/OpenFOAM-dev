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

#include "singleLayerRegion.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "Time.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(singleLayerRegion, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::singleLayerRegion::singleLayerRegion
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, false),
    nHat_
    (
        IOobject
        (
            "nHat",
            time_.timeName(),
            regionMesh()
        ),
        regionMesh(),
        dimensionedVector(dimless, Zero),
        zeroGradientFvPatchField<vector>::typeName
    ),
    magSf_
    (
        IOobject
        (
            "magSf",
            time_.timeName(),
            regionMesh()
        ),
        regionMesh(),
        dimensionedScalar(dimArea, 0)
    ),
    VbyA_
    (
        IOobject
        (
            "VbyA",
            time_.timeName(),
            regionMesh()
        ),
        regionMesh(),
        dimensionedScalar(dimLength, 0),
        zeroGradientFvPatchField<vector>::typeName
    ),
    passivePatchIDs_()
{
    label nBoundaryFaces = 0;
    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& pp = rbm[patchi];
        const labelList& fCells = pp.faceCells();

        nBoundaryFaces += fCells.size();

        UIndirectList<vector>(nHat_, fCells) = pp.faceNormals();
        UIndirectList<scalar>(magSf_, fCells) = pp.magFaceAreas();
    }
    nHat_.correctBoundaryConditions();

    if (nBoundaryFaces != regionMesh().nCells())
    {
        FatalErrorInFunction
            << "Number of primary region coupled boundary faces not equal to "
            << "the number of cells in the local region" << nl << nl
            << "Number of cells = " << regionMesh().nCells() << nl
            << "Boundary faces  = " << nBoundaryFaces << nl
            << abort(FatalError);
    }

    passivePatchIDs_.setSize(intCoupledPatchIDs_.size(), -1);
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppIntCoupled = rbm[patchi];
        if (ppIntCoupled.size() > 0)
        {
            const label cellId = rbm[patchi].faceCells()[0];
            const cell& cFaces = regionMesh().cells()[cellId];
            const label faceO
            (
                cFaces.opposingFaceLabel
                (
                    ppIntCoupled.start(), regionMesh().faces()
                )
            );
            passivePatchIDs_[i] = rbm.whichPatch(faceO);
        }
    }

    Pstream::listCombineGather(passivePatchIDs_, maxEqOp<label>());
    Pstream::listCombineScatter(passivePatchIDs_);

    VbyA_.primitiveFieldRef() = regionMesh().V()/magSf_;
    VbyA_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::singleLayerRegion::~singleLayerRegion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::regionModels::singleLayerRegion::nHat() const
{
    return nHat_;
}


const Foam::volScalarField::Internal&
Foam::regionModels::singleLayerRegion::magSf() const
{
    return magSf_;
}


const Foam::volScalarField& Foam::regionModels::singleLayerRegion::VbyA() const
{
    return VbyA_;
}


const Foam::labelList&
Foam::regionModels::singleLayerRegion::passivePatchIDs() const
{
    return passivePatchIDs_;
}


// ************************************************************************* //
