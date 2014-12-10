/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
        FatalErrorIn("void enrichedPatch::calcMasterPointFaces() const")
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

    Map<DynamicList<label> > mpf(meshPoints().size());

    const faceList& ef = enrichedFaces();

    // Add the original face points
    forAll(masterPatch_, faceI)
    {
        const face& curFace = ef[faceI + slavePatch_.size()];
//         Pout<< "Cur face in pfAddr: " << curFace << endl;
        forAll(curFace, pointI)
        {
            Map<DynamicList<label> >::iterator mpfIter =
                mpf.find(curFace[pointI]);

            if (mpfIter == mpf.end())
            {
                // Not found, add new dynamic list
                mpf.insert
                (
                    curFace[pointI],
                    DynamicList<label>(primitiveMesh::facesPerPoint_)
                );

                // Iterator is invalidated - have to find again
                mpf.find(curFace[pointI])().append(faceI);
            }
            else
            {
                mpfIter().append(faceI);
            }
        }
    }

    // Add the projected points which hit the face
    const labelList& slaveMeshPoints = slavePatch_.meshPoints();

    forAll(slavePointFaceHits_, pointI)
    {
        if
        (
            slavePointPointHits_[pointI] < 0
         && slavePointEdgeHits_[pointI] < 0
         && slavePointFaceHits_[pointI].hit()
        )
        {
            // Get the index of projected point corresponding to this slave
            // point
            const label mergedSmp =
                pointMergeMap().find(slaveMeshPoints[pointI])();

            Map<DynamicList<label> >::iterator mpfIter =
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
                    slavePointFaceHits_[pointI].hitObject()
                );
            }
            else
            {
                mpfIter().append(slavePointFaceHits_[pointI].hitObject());
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
