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

#include "patchProbes.H"
#include "volFields.H"
#include "IOmanip.H"
#include "RemoteData.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchProbes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchProbes,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchProbes::findElements(const fvMesh& mesh)
{
    (void)mesh.tetBasePtIs();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    label patchi = bm.findPatchID(patchName_);

    if (patchi == -1)
    {
        FatalErrorInFunction
            << " Unknown patch name "
            << patchName_ << endl
            << exit(FatalError);
    }

     // All the info for nearest. Construct to miss
    List<RemoteData<scalar>> nearest(this->size());

    const polyPatch& pp = bm[patchi];

    if (pp.size() > 0)
    {
        labelList bndFaces(pp.size());
        forAll(bndFaces, i)
        {
            bndFaces[i] =  pp.start() + i;
        }

        treeBoundBox overallBb(pp.points());
        overallBb = overallBb.extend(1e-4);

        const indexedOctree<treeDataFace> boundaryTree
        (
            treeDataFace    // all information needed to search faces
            (
                false,                      // do not cache bb
                mesh,
                bndFaces                    // patch faces only
            ),
            overallBb,                      // overall search domain
            8,                              // maxLevel
            10,                             // leafsize
            3.0                             // duplicity
        );

        forAll(probeLocations(), probei)
        {
            const point sample = probeLocations()[probei];

            scalar span = boundaryTree.bb().mag();

            pointIndexHit info = boundaryTree.findNearest
            (
                sample,
                Foam::sqr(span)
            );

            if (!info.hit())
            {
                info = boundaryTree.findNearest
                (
                    sample,
                    Foam::sqr(great)
                );
            }

            label facei = boundaryTree.shapes().faceLabels()[info.index()];

            const label patchi = bm.whichPatch(facei);

            if (isA<emptyPolyPatch>(bm[patchi]))
            {
                WarningInFunction
                << " The sample point: " << sample
                << " belongs to " << patchi
                << " which is an empty patch. This is not permitted. "
                << " This sample will not be included "
                << endl;
            }
            else
            {
                const point& fc = mesh.faceCentres()[facei];

                nearest[probei].proci = Pstream::myProcNo();
                nearest[probei].elementi = facei;
                nearest[probei].data = magSqr(fc-sample);
            }
        }
    }

    // Find nearest.
    Pstream::listCombineGather
    (
        nearest,
        RemoteData<scalar>::smallestEqOp()
    );
    Pstream::listCombineScatter(nearest);

    if (debug)
    {
        InfoInFunction << endl;
        forAll(nearest, sampleI)
        {
            Info<< "    " << sampleI << " coord:" << operator[](sampleI)
                << " found on processor:" << nearest[sampleI].proci
                << " in local cell/face:" << nearest[sampleI].elementi
                << endl;
        }
    }

    // Extract any local faces to sample
    elementList_.setSize(nearest.size(), -1);
    forAll(nearest, sampleI)
    {
        if (nearest[sampleI].proci == Pstream::myProcNo())
        {
            elementList_[sampleI] = nearest[sampleI].elementi;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchProbes::patchProbes
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    probes(name, t, dict)
{
    // When constructing probes above it will have called the
    // probes::findElements (since the virtual mechanism not yet operating).
    // Not easy to workaround (apart from feeding through flag into constructor)
    // so clear out any cells found for now.
    elementList_.clear();
    faceList_.clear();

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchProbes::~patchProbes()
{}


bool Foam::patchProbes::write()
{
    if (this->size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }

    return true;
}


bool Foam::patchProbes::read(const dictionary& dict)
{
    dict.lookup("patchName") >> patchName_;
    return probes::read(dict);
}


// ************************************************************************* //
