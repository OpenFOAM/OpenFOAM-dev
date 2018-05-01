/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "regionModel1D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(regionModel1D, 0);
}
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::regionModels::regionModel1D::constructMeshObjects()
{

    nMagSfPtr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "nMagSf",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            regionMesh(),
            dimensionedScalar("zero", dimArea, 0.0)
        )
    );
}


void Foam::regionModels::regionModel1D::initialise()
{
    if (debug)
    {
        Pout<< "regionModel1D::initialise()" << endl;
    }

    // Calculate boundaryFaceFaces and boundaryFaceCells

    DynamicList<label> faceIDs;
    DynamicList<label> cellIDs;

    label localPyrolysisFacei = 0;

    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchi];
        forAll(ppCoupled, localFacei)
        {
            label facei = ppCoupled.start() + localFacei;
            label celli = -1;
            label nFaces = 0;
            label nCells = 0;
            do
            {
                label ownCelli = regionMesh().faceOwner()[facei];
                if (ownCelli != celli)
                {
                    celli = ownCelli;
                }
                else
                {
                    celli = regionMesh().faceNeighbour()[facei];
                }
                nCells++;
                cellIDs.append(celli);
                const cell& cFaces = regionMesh().cells()[celli];
                facei = cFaces.opposingFaceLabel(facei, regionMesh().faces());
                faceIDs.append(facei);
                nFaces++;
            } while (regionMesh().isInternalFace(facei));

            boundaryFaceOppositeFace_[localPyrolysisFacei] = facei;
            faceIDs.remove(); // remove boundary face.
            nFaces--;

            boundaryFaceFaces_[localPyrolysisFacei].transfer(faceIDs);
            boundaryFaceCells_[localPyrolysisFacei].transfer(cellIDs);

            localPyrolysisFacei++;
            nLayers_ = nCells;
        }
    }

    boundaryFaceOppositeFace_.setSize(localPyrolysisFacei);
    boundaryFaceFaces_.setSize(localPyrolysisFacei);
    boundaryFaceCells_.setSize(localPyrolysisFacei);

    surfaceScalarField& nMagSf = nMagSfPtr_();

    surfaceScalarField::Boundary nMagSfBf =
        nMagSf.boundaryFieldRef();

    localPyrolysisFacei = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchi];
        const vectorField& pNormals = ppCoupled.faceNormals();
        nMagSfBf[patchi] = regionMesh().Sf().boundaryField()[patchi] & pNormals;
        forAll(pNormals, localFacei)
        {
            const vector& n = pNormals[localFacei];
            const labelList& faces = boundaryFaceFaces_[localPyrolysisFacei++];
            forAll(faces, facei)
            {
                const label faceID = faces[facei];
                nMagSf[faceID] = regionMesh().Sf()[faceID] & n;
            }
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::regionModels::regionModel1D::read()
{
    if (regionModel::read())
    {
        moveMesh_.readIfPresent("moveMesh", coeffs_);

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::regionModels::regionModel1D::read(const dictionary& dict)
{
    if (regionModel::read(dict))
    {
        moveMesh_.readIfPresent("moveMesh", coeffs_);

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::labelField> Foam::regionModels::regionModel1D::moveMesh
(
    const scalarList& deltaV,
    const scalar minDelta
)
{
    tmp<labelField> tcellMoveMap(new labelField(regionMesh().nCells(), 0));
    labelField& cellMoveMap = tcellMoveMap.ref();

    if (!moveMesh_)
    {
        return cellMoveMap;
    }

    pointField oldPoints = regionMesh().points();
    pointField newPoints = oldPoints;

    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();

    label totalFaceId = 0;
    forAll(intCoupledPatchIDs_, localPatchi)
    {
        label patchi = intCoupledPatchIDs_[localPatchi];
        const polyPatch pp = bm[patchi];
        const vectorField& cf = regionMesh().Cf().boundaryField()[patchi];

        forAll(pp, patchFacei)
        {
            const labelList& faces = boundaryFaceFaces_[totalFaceId];
            const labelList& cells = boundaryFaceCells_[totalFaceId];

            const vector n = pp.faceNormals()[patchFacei];
            const vector sf = pp.faceAreas()[patchFacei];

            List<point> oldCf(faces.size() + 1);
            oldCf[0] = cf[patchFacei];
            forAll(faces, i)
            {
                oldCf[i + 1] = regionMesh().faceCentres()[faces[i]];
            }

            vector newDelta = Zero;
            point nbrCf = oldCf[0];

            forAll(faces, i)
            {
                const label facei = faces[i];
                const label celli = cells[i];

                const face f = regionMesh().faces()[facei];

                newDelta += (deltaV[celli]/mag(sf))*n;

                vector localDelta = Zero;
                forAll(f, pti)
                {
                    const label pointi = f[pti];

                    if
                    (
                        mag((nbrCf - (oldPoints[pointi] + newDelta)) & n)
                      > minDelta
                    )
                    {
                        newPoints[pointi] = oldPoints[pointi] + newDelta;
                        localDelta = newDelta;
                        cellMoveMap[celli] = 1;
                    }
                }
                nbrCf = oldCf[i + 1] + localDelta;
            }
            // Modify boundary
            const label bFacei = boundaryFaceOppositeFace_[totalFaceId];
            const face f = regionMesh().faces()[bFacei];
            const label celli = cells[cells.size() - 1];
            newDelta += (deltaV[celli]/mag(sf))*n;
            forAll(f, pti)
            {
                const label pointi = f[pti];
                if
                (
                    mag((nbrCf - (oldPoints[pointi] + newDelta)) & n)
                  > minDelta
                )
                {
                    newPoints[pointi] = oldPoints[pointi] + newDelta;
                    cellMoveMap[celli] = 1;
                }
            }
            totalFaceId ++;
        }
    }
    // Move points
    regionMesh().movePoints(newPoints);

    return tcellMoveMap;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType
)
:
    regionModel(mesh, regionType),
    boundaryFaceFaces_(),
    boundaryFaceCells_(),
    boundaryFaceOppositeFace_(),
    nLayers_(0),
    nMagSfPtr_(nullptr),
    moveMesh_(false)
{}


Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, false),
    boundaryFaceFaces_(regionMesh().nCells()),
    boundaryFaceCells_(regionMesh().nCells()),
    boundaryFaceOppositeFace_(regionMesh().nCells()),
    nLayers_(0),
    nMagSfPtr_(nullptr),
    moveMesh_(true)
{
    if (active_)
    {
        constructMeshObjects();
        initialise();

        if (readFields)
        {
            read();
        }
    }
}


Foam::regionModels::regionModel1D::regionModel1D
(
    const fvMesh& mesh,
    const word& regionType,
    const word& modelName,
    const dictionary& dict,
    bool readFields
)
:
    regionModel(mesh, regionType, modelName, dict, readFields),
    boundaryFaceFaces_(regionMesh().nCells()),
    boundaryFaceCells_(regionMesh().nCells()),
    boundaryFaceOppositeFace_(regionMesh().nCells()),
    nLayers_(0),
    nMagSfPtr_(nullptr),
    moveMesh_(false)
{
    if (active_)
    {
        constructMeshObjects();
        initialise();

        if (readFields)
        {
            read(dict);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::regionModel1D::~regionModel1D()
{}


// ************************************************************************* //
