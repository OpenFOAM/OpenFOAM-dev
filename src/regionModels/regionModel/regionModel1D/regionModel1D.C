/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

    label localPyrolysisFaceI = 0;

    const polyBoundaryMesh& rbm = regionMesh().boundaryMesh();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchI];
        forAll(ppCoupled, localFaceI)
        {
            label faceI = ppCoupled.start() + localFaceI;
            label cellI = -1;
            label nFaces = 0;
            label nCells = 0;
            do
            {
                label ownCellI = regionMesh().faceOwner()[faceI];
                if (ownCellI != cellI)
                {
                    cellI = ownCellI;
                }
                else
                {
                    cellI = regionMesh().faceNeighbour()[faceI];
                }
                nCells++;
                cellIDs.append(cellI);
                const cell& cFaces = regionMesh().cells()[cellI];
                faceI = cFaces.opposingFaceLabel(faceI, regionMesh().faces());
                faceIDs.append(faceI);
                nFaces++;
            } while (regionMesh().isInternalFace(faceI));

            boundaryFaceOppositeFace_[localPyrolysisFaceI] = faceI;
            faceIDs.remove(); //remove boundary face.
            nFaces--;

            boundaryFaceFaces_[localPyrolysisFaceI].transfer(faceIDs);
            boundaryFaceCells_[localPyrolysisFaceI].transfer(cellIDs);

            localPyrolysisFaceI++;
            nLayers_ = nCells;
        }
    }

    boundaryFaceOppositeFace_.setSize(localPyrolysisFaceI);
    boundaryFaceFaces_.setSize(localPyrolysisFaceI);
    boundaryFaceCells_.setSize(localPyrolysisFaceI);

    surfaceScalarField& nMagSf = nMagSfPtr_();

    localPyrolysisFaceI = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchI = intCoupledPatchIDs_[i];
        const polyPatch& ppCoupled = rbm[patchI];
        const vectorField& pNormals = ppCoupled.faceNormals();
        nMagSf.boundaryField()[patchI] =
            regionMesh().Sf().boundaryField()[patchI] & pNormals;
        forAll(pNormals, localFaceI)
        {
            const vector& n = pNormals[localFaceI];
            const labelList& faces = boundaryFaceFaces_[localPyrolysisFaceI++];
            forAll (faces, faceI)
            {
                const label faceID = faces[faceI];
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
    labelField& cellMoveMap = tcellMoveMap();

    if (!moveMesh_)
    {
        return cellMoveMap;
    }

    pointField oldPoints = regionMesh().points();
    pointField newPoints = oldPoints;

    const polyBoundaryMesh& bm = regionMesh().boundaryMesh();

    label totalFaceId = 0;
    forAll(intCoupledPatchIDs_, localPatchI)
    {
        label patchI = intCoupledPatchIDs_[localPatchI];
        const polyPatch pp = bm[patchI];
        const vectorField& cf = regionMesh().Cf().boundaryField()[patchI];

        forAll(pp, patchFaceI)
        {
            const labelList& faces = boundaryFaceFaces_[totalFaceId];
            const labelList& cells = boundaryFaceCells_[totalFaceId];

            const vector n = pp.faceNormals()[patchFaceI];
            const vector sf = pp.faceAreas()[patchFaceI];

            List<point> oldCf(faces.size() + 1);
            oldCf[0] = cf[patchFaceI];
            forAll(faces, i)
            {
                oldCf[i + 1] = regionMesh().faceCentres()[faces[i]];
            }

            vector newDelta = vector::zero;
            point nbrCf = oldCf[0];

            forAll(faces, i)
            {
                const label faceI = faces[i];
                const label cellI = cells[i];

                const face f = regionMesh().faces()[faceI];

                newDelta += (deltaV[cellI]/mag(sf))*n;

                vector localDelta = vector::zero;
                forAll(f, pti)
                {
                    const label pointI = f[pti];

                    if
                    (
                        mag((nbrCf - (oldPoints[pointI] + newDelta)) & n)
                      > minDelta
                    )
                    {
                        newPoints[pointI] = oldPoints[pointI] + newDelta;
                        localDelta = newDelta;
                        cellMoveMap[cellI] = 1;
                    }
                }
                nbrCf = oldCf[i + 1] + localDelta;
            }
            // Modify boundary
            const label bFaceI = boundaryFaceOppositeFace_[totalFaceId];
            const face f = regionMesh().faces()[bFaceI];
            const label cellI = cells[cells.size() - 1];
            newDelta += (deltaV[cellI]/mag(sf))*n;
            forAll(f, pti)
            {
                const label pointI = f[pti];
                if
                (
                    mag((nbrCf - (oldPoints[pointI] + newDelta)) & n)
                  > minDelta
                )
                {
                    newPoints[pointI] = oldPoints[pointI] + newDelta;
                    cellMoveMap[cellI] = 1;
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
    nMagSfPtr_(NULL),
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
    nMagSfPtr_(NULL),
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
    nMagSfPtr_(NULL),
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
