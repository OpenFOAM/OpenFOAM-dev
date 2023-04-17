/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "patchInjectionBase.H"
#include "polyMesh.H"
#include "SubField.H"
#include "Random.H"
#include "triPointRef.H"
#include "volFields.H"
#include "polyMeshTetDecomposition.H"
#include "polygonTriangulate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchInjectionBase::patchInjectionBase
(
    const polyMesh& mesh,
    const word& patchName
)
:
    patchName_(patchName),
    patchId_(mesh.boundaryMesh().findPatchID(patchName_)),
    sumProcArea_(),
    sumFaceArea_(),
    sumFaceTriArea_()
{
    if (patchId_ < 0)
    {
        FatalErrorInFunction
            << "Requested patch " << patchName_ << " not found" << nl
            << "Available patches are: " << mesh.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    topoChange(mesh);
}


Foam::patchInjectionBase::patchInjectionBase(const patchInjectionBase& pib)
:
    patchName_(pib.patchName_),
    patchId_(pib.patchId_),
    sumProcArea_(pib.sumProcArea_),
    sumFaceArea_(pib.sumFaceArea_),
    sumFaceTriArea_(pib.sumFaceTriArea_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchInjectionBase::~patchInjectionBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchInjectionBase::topoChange(const polyMesh& mesh)
{
    // Set/cache the injector cells
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    // Initialise
    sumProcArea_.resize(Pstream::nProcs());
    sumProcArea_ = 0;
    sumFaceArea_.resize(patch.size());
    sumFaceArea_ = 0;
    sumFaceTriArea_.resize(patch.size());
    forAll(patch, patchFacei)
    {
        sumFaceTriArea_[patchFacei].resize(patch[patchFacei].nTriangles());
        sumFaceTriArea_[patchFacei] = 0;
    }

    // Cumulatively sum up the areas through the local patch
    scalar patchArea = 0;
    forAll(patch, patchFacei)
    {
        const label facei = patchFacei + patch.start();
        const label celli = patch.faceCells()[patchFacei];

        scalar patchFaceArea = 0;
        for
        (
            label patchFaceTrii = 0;
            patchFaceTrii < patch[patchFacei].nTriangles();
            ++ patchFaceTrii
        )
        {
            const tetIndices tet(celli, facei, patchFaceTrii + 1);

            patchFaceArea += tet.faceTri(mesh).mag();
            sumFaceTriArea_[patchFacei][patchFaceTrii] = patchFaceArea;
        }

        patchArea += patchFaceArea;
        sumFaceArea_[patchFacei] = patchArea;
    }

    // Cumulatively sum the total areas across the processors
    sumProcArea_[Pstream::myProcNo()] = patchArea;
    Pstream::listCombineGather(sumProcArea_, maxEqOp<scalar>());
    Pstream::listCombineScatter(sumProcArea_);
    for (label proci = 1; proci < Pstream::nProcs(); proci ++)
    {
        sumProcArea_[proci] += sumProcArea_[proci - 1];
    }
}


void Foam::patchInjectionBase::setPositionAndCell
(
    const fvMesh& mesh,
    Random& rnd,
    barycentric& coordinates,
    label& celli,
    label& tetFacei,
    label& tetPti,
    label& facei
)
{
    // Linear searching for area fractions
    auto findArea = []
    (
        const scalarList& areas,
        scalar& area
    )
    {
        for (label i = areas.size(); i > 0; i --)
        {
            if (area >= areas[i - 1])
            {
                area -= areas[i - 1];
                return i;
            }
        }

        return label(0);
    };

    const polyPatch& patch = mesh.boundaryMesh()[patchId_];

    scalar area = rnd.globalScalar01()*sumProcArea_.last();

    if (patch.size() > 0)
    {
        // Determine which processor to inject from
        const label proci = findArea(sumProcArea_, area);

        // If this processor...
        if (Pstream::myProcNo() == proci)
        {
            // Determine the face to inject from
            const label patchFacei = findArea(sumFaceArea_, area);

            // Determine the face-tri to inject from
            const label patchFaceTrii =
                findArea(sumFaceTriArea_[patchFacei], area);

            // Set the topology
            const barycentric2D r = barycentric2D01(rnd);
            coordinates = barycentric(0, r.a(), r.b(), r.c());
            celli = mesh.faceOwner()[patch.start() + patchFacei];
            tetFacei = patch.start() + patchFacei;
            tetPti = patchFaceTrii + 1;
            facei = patch.start() + patchFacei;
        }
        else
        {
            coordinates = barycentric::uniform(NaN);
            celli = -1;
            tetFacei = -1;
            tetPti = -1;
            facei = -1;
        }
    }
    else
    {
        coordinates = barycentric::uniform(NaN);
        celli = -1;
        tetFacei = -1;
        tetPti = -1;
        facei = -1;
    }
}


// ************************************************************************* //
