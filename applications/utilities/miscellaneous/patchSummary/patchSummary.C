/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Application
    patchSummary

Description
    Writes fields and boundary condition info for each patch at each requested
    time instance.

    Default action is to write a single entry for patches/patchGroups with the
    same boundary conditions. Use the -expand option to print every patch
    separately. In case of multiple groups matching it will print only the
    first one.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "pointFields.H"
#include "IOobjectList.H"
#include "patchSummaryTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "addRegionOption.H"
    argList::addBoolOption
    (
        "expand",
        "Do not combine patches"
    );
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    const bool expand = args.optionFound("expand");


#   include "createNamedMesh.H"
    const polyBoundaryMesh& bm = mesh.boundaryMesh();


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Update the mesh if changed
        if (mesh.readUpdate() == polyMesh::TOPO_PATCH_CHANGE)
        {
            Info<< "Detected changed patches. Recreating patch group table."
                << endl;
        }


        const IOobjectList fieldObjs(mesh, runTime.timeName());
        const wordList objNames = fieldObjs.names();

        PtrList<volScalarField> vsf(objNames.size());
        PtrList<volVectorField> vvf(objNames.size());
        PtrList<volSphericalTensorField> vsptf(objNames.size());
        PtrList<volSymmTensorField> vsytf(objNames.size());
        PtrList<volTensorField> vtf(objNames.size());

        PtrList<pointScalarField> psf(objNames.size());
        PtrList<pointVectorField> pvf(objNames.size());
        PtrList<pointSphericalTensorField> psptf(objNames.size());
        PtrList<pointSymmTensorField> psytf(objNames.size());
        PtrList<pointTensorField> ptf(objNames.size());

        Info<< "Valid fields:" << endl;

        forAll(objNames, objI)
        {
            IOobject obj
            (
                objNames[objI],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            if (obj.headerOk())
            {
                addToFieldList(vsf, obj, objI, mesh);
                addToFieldList(vvf, obj, objI, mesh);
                addToFieldList(vsptf, obj, objI, mesh);
                addToFieldList(vsytf, obj, objI, mesh);
                addToFieldList(vtf, obj, objI, mesh);

                addToFieldList(psf, obj, objI, pointMesh::New(mesh));
                addToFieldList(pvf, obj, objI, pointMesh::New(mesh));
                addToFieldList(psptf, obj, objI, pointMesh::New(mesh));
                addToFieldList(psytf, obj, objI, pointMesh::New(mesh));
                addToFieldList(ptf, obj, objI, pointMesh::New(mesh));
            }
        }

        Info<< endl;


        if (expand)
        {
            // Print each patch separately

            forAll(bm, patchI)
            {
                Info<< bm[patchI].type() << "\t: " << bm[patchI].name() << nl;
                outputFieldList(vsf, patchI);
                outputFieldList(vvf, patchI);
                outputFieldList(vsptf, patchI);
                outputFieldList(vsytf, patchI);
                outputFieldList(vtf, patchI);

                outputFieldList(psf, patchI);
                outputFieldList(pvf, patchI);
                outputFieldList(psptf, patchI);
                outputFieldList(psytf, patchI);
                outputFieldList(ptf, patchI);
                Info<< endl;
            }
        }
        else
        {
            // Collect for each patch the bc type per field. Merge similar
            // patches.

            // Per 'group', the map from fieldname to patchfield type
            DynamicList<HashTable<word> > fieldToTypes(bm.size());
            // Per 'group' the patches
            DynamicList<DynamicList<label> > groupToPatches(bm.size());
            forAll(bm, patchI)
            {
                HashTable<word> fieldToType;
                collectFieldList(vsf, patchI, fieldToType);
                collectFieldList(vvf, patchI, fieldToType);
                collectFieldList(vsptf, patchI, fieldToType);
                collectFieldList(vsytf, patchI, fieldToType);
                collectFieldList(vtf, patchI, fieldToType);

                collectFieldList(psf, patchI, fieldToType);
                collectFieldList(pvf, patchI, fieldToType);
                collectFieldList(psptf, patchI, fieldToType);
                collectFieldList(psytf, patchI, fieldToType);
                collectFieldList(ptf, patchI, fieldToType);

                label groupI = findIndex(fieldToTypes, fieldToType);
                if (groupI == -1)
                {
                    DynamicList<label> group(1);
                    group.append(patchI);
                    groupToPatches.append(group);
                    fieldToTypes.append(fieldToType);
                }
                else
                {
                    groupToPatches[groupI].append(patchI);
                }
            }


            forAll(groupToPatches, groupI)
            {
                const DynamicList<label>& patchIDs = groupToPatches[groupI];

                if (patchIDs.size() > 1)
                {
                    // Check if part of a group
                    wordList groups;
                    labelHashSet nonGroupPatches;
                    bm.matchGroups(patchIDs, groups, nonGroupPatches);

                    const labelList sortedPatches(nonGroupPatches.sortedToc());
                    forAll(sortedPatches, i)
                    {
                        Info<< bm[sortedPatches[i]].type()
                            << "\t: " << bm[sortedPatches[i]].name() << nl;
                    }
                    if (groups.size())
                    {
                        forAll(groups, i)
                        {
                            Info<< "group\t: " << groups[i] << nl;
                        }
                    }
                    outputFieldList(vsf, patchIDs[0]);
                    outputFieldList(vvf, patchIDs[0]);
                    outputFieldList(vsptf, patchIDs[0]);
                    outputFieldList(vsytf, patchIDs[0]);
                    outputFieldList(vtf, patchIDs[0]);

                    outputFieldList(psf, patchIDs[0]);
                    outputFieldList(pvf, patchIDs[0]);
                    outputFieldList(psptf, patchIDs[0]);
                    outputFieldList(psytf, patchIDs[0]);
                    outputFieldList(ptf, patchIDs[0]);
                    Info<< endl;
                }
                else
                {
                    // No group.
                    forAll(patchIDs, i)
                    {
                        label patchI = patchIDs[i];
                        Info<< bm[patchI].type()
                            << "\t: " << bm[patchI].name() << nl;
                        outputFieldList(vsf, patchI);
                        outputFieldList(vvf, patchI);
                        outputFieldList(vsptf, patchI);
                        outputFieldList(vsytf, patchI);
                        outputFieldList(vtf, patchI);

                        outputFieldList(psf, patchI);
                        outputFieldList(pvf, patchI);
                        outputFieldList(psptf, patchI);
                        outputFieldList(psytf, patchI);
                        outputFieldList(ptf, patchI);
                        Info<< endl;
                    }
                }
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
