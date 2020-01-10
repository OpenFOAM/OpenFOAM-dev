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

Description
    Create intermediate mesh files from SAMM files

\*---------------------------------------------------------------------------*/

#include "sammMesh.H"
#include "Time.H"
#include "wallPolyPatch.H"
#include "mergedCyclicPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "preservePatchTypes.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::sammMesh::readBoundary()
{
    label nPatches=0, nFaces=0;
    labelList nPatchFaces(1000);

    label lineIndex, sammLabel;
    label sammRegion, configNumber;

    labelList pointLabels(4);
    labelList pointLabelsTri(3);

    labelList patchLabels(1000, label(-1));

    word patchType;
    patchTypes_.setSize(1000);
    patchNames_.setSize(1000);

    fileName boundaryFileName(casePrefix_ + ".bnd");

    {
        IFstream boundaryFile(boundaryFileName);

        // Collect no. of faces (nFaces),
        // no. of patches (nPatches)
        // and for each of these patches the number of faces
        // (nPatchFaces[patchLabel])
        // and a conversion table from Samm regions to (Foam) patchLabels

        if (boundaryFile.good())
        {
            forAll(nPatchFaces, faceLabel)
            {
                nPatchFaces[faceLabel] = 0;
            }

            while ((boundaryFile >> lineIndex).good())
            {
                nFaces++;

                // Skip point numbers
                for (int i=0; i<4; i++)
                {
                    boundaryFile >> sammLabel;
                }

                boundaryFile >> sammRegion;
                boundaryFile >> configNumber;
                boundaryFile >> patchType;

                // Build translation table to convert samm patch to foam patch
                label patchLabel = patchLabels[sammRegion];
                if (patchLabel == -1)
                {
                    patchLabel = nPatches;
                    patchLabels[sammRegion] = patchLabel;
                    patchTypes_[patchLabel] = patchType;
                    patchNames_[patchLabel] = patchType + name(sammRegion);

                    nPatches++;

                    Info<< "Samm region " << sammRegion
                        << " with type " << patchType
                        << " is now Foam patch " << patchLabel << endl;

                }

                nPatchFaces[patchLabel]++;
            }


            Info<< nl
                << "Setting size of shapePatchList to " << nPatches
                << nl << endl;

            nPatchFaces.setSize(nPatches);
            patchTypes_.setSize(nPatches);
            patchNames_.setSize(nPatches);
        }
        else
        {
            FatalErrorInFunction
                << "Cannot read file "
                << boundaryFileName
                << abort(FatalError);
        }
    }

    if (nPatches > 0)
    {
        boundary_.setSize(nPatchFaces.size());
        patchTypes_.setSize(nPatchFaces.size());
        patchNames_.setSize(nPatchFaces.size());

        // size the lists and reset the counters to be used again
        forAll(boundary_, patchLabel)
        {
            boundary_[patchLabel].setSize(nPatchFaces[patchLabel]);

            nPatchFaces[patchLabel] = 0;
        }

        IFstream boundaryFile(boundaryFileName);

        for (label facei=0; facei<nFaces; facei++)
        {
            boundaryFile >> lineIndex;

            for (int i = 0; i < 4; i++)
            {
                boundaryFile >> sammLabel;

                // convert Samm label to Foam point label
                // through lookup-list starPointLabelLookup_
                pointLabels[i] = starPointLabelLookup_[sammLabel];

                if (pointLabels[i] < 0)
                {
                    Info<< "Boundary file not consistent with vertex file\n"
                        << "Samm vertex number " << sammLabel
                        << " does not exist\n";
                }

            }

            boundaryFile >> sammRegion;
            label patchLabel = patchLabels[sammRegion];

            boundaryFile >> configNumber;
            boundaryFile >> patchType;

            if   // Triangle
            (
                pointLabels[2] == pointLabels[3]
            )
            {
                // Info<< "Converting collapsed quad into triangle"
                //    << " for face " << facei
                //    << " in Samm boundary " << lineIndex << endl;

                pointLabelsTri[0] = pointLabels[0];
                pointLabelsTri[1] = pointLabels[1];
                pointLabelsTri[2] = pointLabels[2];

                boundary_[patchLabel][nPatchFaces[patchLabel]]
                    = face(pointLabelsTri);
            }
            else
            {
                boundary_[patchLabel][nPatchFaces[patchLabel]]
                    = face(pointLabels);
            }

            // increment counter of faces in current patch
            nPatchFaces[patchLabel]++;
        }

        forAll(boundary_, patchLabel)
        {
            word patchType = patchTypes_[patchLabel];

            if (patchType == "SYMP")
            {
                patchTypes_[patchLabel] = symmetryPolyPatch::typeName;
            }
            else if (patchType == "WALL")
            {
                patchTypes_[patchLabel] = wallPolyPatch::typeName;
            }
            else if (patchType == "CYCL")
            {
                patchTypes_[patchLabel] = mergedCyclicPolyPatch::typeName;
            }
            else
            {
                patchTypes_[patchLabel] = polyPatch::typeName;
            }

            Info<< "Foam patch " << patchLabel
                << " is of type " << patchTypes_[patchLabel]
                << " with name " << patchNames_[patchLabel] << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "No boundary faces in file "
            << boundaryFileName
            << endl;
    }

    patchPhysicalTypes_.setSize(patchTypes_.size());

    PtrList<dictionary> patchDicts;

    preservePatchTypes
    (
        runTime_,
        runTime_.constant(),
        polyMesh::meshSubDir,
        patchNames_,
        patchDicts,
        defaultFacesName_,
        defaultFacesType_
    );

    forAll(patchDicts, patchi)
    {
        if (patchDicts.set(patchi))
        {
            const dictionary& dict = patchDicts[patchi];
            dict.readIfPresent("type", patchTypes_[patchi]);
            dict.readIfPresent("physicalType", patchPhysicalTypes_[patchi]);
        }
    }
}


// ************************************************************************* //
