/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "enrichedPatch.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::enrichedPatch::nFaceHits_ = 4;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcMasterPointFaces() const
{
    if (masterPointFacesPtr_)
    {
        FatalErrorInFunction
            << "Master point face addressing already calculated."
            << abort(FatalError);
    }

    // Note:
    // Master point face addressing lists the master faces for all points
    // in the enriched patch support (if there are no master faces, which is
    // normal, the list will be empty).  The index represents the index of
    // the master face rather than the index from the enriched patch
    // Master face points lists the points of the enriched master face plus
    // points projected into the master face

    Map<DynamicList<label>> mpf(meshPoints().size());

    const faceList& ef = enrichedFaces();

    // Add the original face points
    forAll(masterPatch_, facei)
    {
        const face& curFace = ef[facei + slavePatch_.size()];
//         Pout<< "Cur face in pfAddr: " << curFace << endl;
        forAll(curFace, pointi)
        {
            Map<DynamicList<label>>::iterator mpfIter =
                mpf.find(curFace[pointi]);

            if (mpfIter == mpf.end())
            {
                // Not found, add new dynamic list
                mpf.insert
                (
                    curFace[pointi],
                    DynamicList<label>(primitiveMesh::facesPerPoint_)
                );

                // Iterator is invalidated - have to find again
                mpf.find(curFace[pointi])().append(facei);
            }
            else
            {
                mpfIter().append(facei);
            }
        }
    }

    // Add the projected points which hit the face
    const labelList& slaveMeshPoints = slavePatch_.meshPoints();

    forAll(slavePointFaceHits_, pointi)
    {
        if
        (
            slavePointPointHits_[pointi] < 0
         && slavePointEdgeHits_[pointi] < 0
         && slavePointFaceHits_[pointi].hit()
        )
        {
            // Get the index of projected point corresponding to this slave
            // point
            const label mergedSmp =
                pointMergeMap().find(slaveMeshPoints[pointi])();

            Map<DynamicList<label>>::iterator mpfIter =
                mpf.find(mergedSmp);

            if (mpfIter == mpf.end())
            {
                // Not found, add new dynamic list
                mpf.insert
                (
                    mergedSmp,
                    DynamicList<label>(primitiveMesh::facesPerPoint_)
                );

                // Iterator is invalidated - have to find again
                mpf.find(mergedSmp)().append
                (
                    slavePointFaceHits_[pointi].hitObject()
                );
            }
            else
            {
                mpfIter().append(slavePointFaceHits_[pointi].hitObject());
            }
        }
    }

    // Re-pack dynamic lists into normal lists
    const labelList mpfToc = mpf.toc();

    masterPointFacesPtr_ = new Map<labelList>(2*mpfToc.size());
    Map<labelList>& masterPointFaceAddr = *masterPointFacesPtr_;

    forAll(mpfToc, mpfTocI)
    {
        labelList l;
        l.transfer(mpf.find(mpfToc[mpfTocI])());

        masterPointFaceAddr.insert(mpfToc[mpfTocI], l);
    }
    // Pout<< "masterPointFaceAddr: " << masterPointFaceAddr << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Map<Foam::labelList>& Foam::enrichedPatch::masterPointFaces() const
{
    if (!masterPointFacesPtr_)
    {
        calcMasterPointFaces();
    }

    return *masterPointFacesPtr_;
}


// ************************************************************************* //
