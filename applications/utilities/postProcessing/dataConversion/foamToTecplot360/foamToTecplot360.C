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

Application
    foamToTecplot360

Description
    Tecplot binary file format writer.

Usage
    \b foamToTecplot360 [OPTION]

    Options:
      - \par -fields \<names\>
        Convert selected fields only. For example,
        \verbatim
          -fields '( p T U )'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument.

      - \par -cellSet \<name\>

      - \par -faceSet \<name\>
        Restrict conversion to the cellSet, faceSet.

      - \par -nearCellValue
        Output cell value on patches instead of patch value itself

      - \par -noInternal
        Do not generate file for mesh, only for patches

      - \par -noPointValues
        No pointFields

      - \par -noFaceZones
        No faceZones

      - \par -excludePatches \<patchNames\>
        Specify patches (wildcards) to exclude. For example,
        \verbatim
          -excludePatches '( inlet_1 inlet_2 "proc.*")'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument. The double quotes denote a regular
        expression.

      - \par -useTimeName
        use the time index in the VTK file name instead of the time index

\*---------------------------------------------------------------------------*/

#include "pointMesh.H"
#include "volPointInterpolation.H"
#include "emptyPolyPatch.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "sphericalTensorIOField.H"
#include "symmTensorIOField.H"
#include "tensorIOField.H"
#include "passiveParticleCloud.H"
#include "faceSet.H"
#include "stringListOps.H"
#include "wordRe.H"

#include "vtkMesh.H"
#include "readFields.H"
#include "tecplotWriter.H"

#include "TECIO.h"

// Note: needs to be after TECIO to prevent Foam::Time conflicting with
// Xlib Time.
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void print(const char* msg, Ostream& os, const PtrList<GeoField>& flds)
{
    if (flds.size())
    {
        os  << msg;
        forAll(flds, i)
        {
            os  << ' ' << flds[i].name();
        }
        os  << endl;
    }
}


void print(Ostream& os, const wordList& flds)
{
    forAll(flds, i)
    {
        os  << ' ' << flds[i];
    }
    os  << endl;
}


labelList getSelectedPatches
(
    const polyBoundaryMesh& patches,
    const List<wordRe>& excludePatches  // HashSet<word>& excludePatches
)
{
    DynamicList<label> patchIDs(patches.size());

    Info<< "Combining patches:" << endl;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
            isType<emptyPolyPatch>(pp)
            || (Pstream::parRun() && isType<processorPolyPatch>(pp))
        )
        {
            Info<< "    discarding empty/processor patch " << patchi
                << " " << pp.name() << endl;
        }
        else if (findStrings(excludePatches, pp.name()))
        {
            Info<< "    excluding patch " << patchi
                << " " << pp.name() << endl;
        }
        else
        {
            patchIDs.append(patchi);
            Info<< "    patch " << patchi << " " << pp.name() << endl;
        }
    }
    return patchIDs.shrink();
}





int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Tecplot binary file format writer"
    );

    timeSelector::addOptions();
    #include "addRegionOption.H"

    argList::addOption
    (
        "fields",
        "names",
        "convert selected fields only. eg, '(p T U)'"
    );
    argList::addOption
    (
        "cellSet",
        "name",
        "restrict conversion to the specified cellSet"
    );
    argList::addOption
    (
        "faceSet",
        "name",
        "restrict conversion to the specified cellSet"
    );
    argList::addBoolOption
    (
        "nearCellValue",
        "output cell value on patches instead of patch value itself"
    );
    argList::addBoolOption
    (
        "noInternal",
        "do not generate file for mesh, only for patches"
    );
    argList::addBoolOption
    (
        "noPointValues",
        "no pointFields"
    );
    argList::addOption
    (
        "excludePatches",
        "patches (wildcards) to exclude"
    );
    argList::addBoolOption
    (
        "noFaceZones",
        "no faceZones"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const bool doWriteInternal = !args.optionFound("noInternal");
    const bool doFaceZones     = !args.optionFound("noFaceZones");
    const bool nearCellValue = args.optionFound("nearCellValue");
    const bool noPointValues = args.optionFound("noPointValues");

    if (nearCellValue)
    {
        WarningInFunction
            << "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }

    if (noPointValues)
    {
        WarningInFunction
            << "Outputting cell values only" << nl << endl;
    }

    List<wordRe> excludePatches;
    if (args.optionFound("excludePatches"))
    {
        args.optionLookup("excludePatches")() >> excludePatches;

        Info<< "Not including patches " << excludePatches << nl << endl;
    }

    word cellSetName;
    string vtkName;

    if (args.optionReadIfPresent("cellSet", cellSetName))
    {
        vtkName = cellSetName;
    }
    else if (Pstream::parRun())
    {
        // Strip off leading casename, leaving just processor_DDD ending.
        vtkName = runTime.caseName();

        string::size_type i = vtkName.rfind("processor");

        if (i != string::npos)
        {
            vtkName = vtkName.substr(i);
        }
    }
    else
    {
        vtkName = runTime.caseName();
    }


    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    // TecplotData/ directory in the case
    fileName fvPath(runTime.path()/"Tecplot360");
    // Directory of mesh (region0 gets filtered out)
    fileName regionPrefix = "";

    if (regionName != polyMesh::defaultRegion)
    {
        fvPath = fvPath/regionName;
        regionPrefix = regionName;
    }

    if (isDir(fvPath))
    {
        if
        (
            args.optionFound("time")
         || args.optionFound("latestTime")
         || cellSetName.size()
         || regionName != polyMesh::defaultRegion
        )
        {
            Info<< "Keeping old files in " << fvPath << nl << endl;
        }
        else
        {
            Info<< "Deleting old VTK files in " << fvPath << nl << endl;

            rmDir(fvPath);
        }
    }

    mkDir(fvPath);


    // mesh wrapper; does subsetting and decomposition
    vtkMesh vMesh(mesh, cellSetName);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time: " << runTime.timeName() << endl;

        const word timeDesc = name(timeI);    // name(runTime.timeIndex());

        // Check for new polyMesh/ and update mesh, fvMeshSubset and cell
        // decomposition.
        polyMesh::readUpdateState meshState = vMesh.readUpdate();

        const fvMesh& mesh = vMesh.mesh();

        INTEGER4 nFaceNodes = 0;
        forAll(mesh.faces(), facei)
        {
            nFaceNodes += mesh.faces()[facei].size();
        }


        // Read all fields on the new mesh
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());

        HashSet<word> selectedFields;
        if (args.optionFound("fields"))
        {
            args.optionLookup("fields")() >> selectedFields;
        }

        // Construct the vol fields (on the original mesh if subsetted)

        PtrList<volScalarField> vsf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vsf);
        print("    volScalarFields            :", Info, vsf);

        PtrList<volVectorField> vvf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vvf);
        print("    volVectorFields            :", Info, vvf);

        PtrList<volSphericalTensorField> vSpheretf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vSpheretf);
        print("    volSphericalTensorFields   :", Info, vSpheretf);

        PtrList<volSymmTensorField> vSymmtf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vSymmtf);
        print("    volSymmTensorFields        :", Info, vSymmtf);

        PtrList<volTensorField> vtf;
        readFields(vMesh, vMesh.baseMesh(), objects, selectedFields, vtf);
        print("    volTensorFields            :", Info, vtf);


        // Construct pointMesh only if necessary since constructs edge
        // addressing (expensive on polyhedral meshes)
        if (noPointValues)
        {
            Info<< "    pointScalarFields : switched off"
                << " (\"-noPointValues\" (at your option)\n";
            Info<< "    pointVectorFields : switched off"
                << " (\"-noPointValues\" (at your option)\n";
        }

        PtrList<pointScalarField> psf;
        PtrList<pointVectorField> pvf;
        // PtrList<pointSphericalTensorField> pSpheretf;
        // PtrList<pointSymmTensorField> pSymmtf;
        // PtrList<pointTensorField> ptf;


        if (!noPointValues)
        {
            //// Add interpolated volFields
            // const volPointInterpolation& pInterp = volPointInterpolation::New
            //(
            //    mesh
            //);
            //
            // label nPsf = psf.size();
            // psf.setSize(nPsf+vsf.size());
            // forAll(vsf, i)
            //{
            //    Info<< "Interpolating " << vsf[i].name() << endl;
            //    tmp<pointScalarField> tvsf(pInterp.interpolate(vsf[i]));
            //    tvsf().rename(vsf[i].name() + "_point");
            //    psf.set(nPsf, tvsf);
            //    nPsf++;
            //}
            //
            // label nPvf = pvf.size();
            // pvf.setSize(vvf.size());
            // forAll(vvf, i)
            //{
            //    Info<< "Interpolating " << vvf[i].name() << endl;
            //    tmp<pointVectorField> tvvf(pInterp.interpolate(vvf[i]));
            //    tvvf().rename(vvf[i].name() + "_point");
            //    pvf.set(nPvf, tvvf);
            //    nPvf++;
            //}

            readFields
            (
                vMesh,
                pointMesh::New(vMesh.baseMesh()),
                objects,
                selectedFields,
                psf
            );
            print("    pointScalarFields          :", Info, psf);

            readFields
            (
                vMesh,
                pointMesh::New(vMesh.baseMesh()),
                objects,
                selectedFields,
                pvf
            );
            print("    pointVectorFields          :", Info, pvf);

            // readFields
            //(
            //    vMesh,
            //    pointMesh::New(vMesh.baseMesh()),
            //    objects,
            //    selectedFields,
            //    pSpheretf
            //);
            // print("    pointSphericalTensorFields :", Info, pSpheretf);
            //
            // readFields
            //(
            //    vMesh,
            //    pointMesh::New(vMesh.baseMesh()),
            //    objects,
            //    selectedFields,
            //    pSymmtf
            //);
            // print("    pointSymmTensorFields      :", Info, pSymmtf);
            //
            // readFields
            //(
            //    vMesh,
            //    pointMesh::New(vMesh.baseMesh()),
            //    objects,
            //    selectedFields,
            //    ptf
            //);
            // print("    pointTensorFields          :", Info, ptf);
        }
        Info<< endl;


        // Get field names
        // ~~~~~~~~~~~~~~~

        string varNames;
        DynamicList<INTEGER4> varLocation;

        string cellVarNames;
        DynamicList<INTEGER4> cellVarLocation;

        // volFields
        tecplotWriter::getTecplotNames
        (
            vsf,
            ValueLocation_CellCentered,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vsf,
            ValueLocation_CellCentered,
            cellVarNames,
            cellVarLocation
        );

        tecplotWriter::getTecplotNames
        (
            vvf,
            ValueLocation_CellCentered,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vvf,
            ValueLocation_CellCentered,
            cellVarNames,
            cellVarLocation
        );

        tecplotWriter::getTecplotNames
        (
            vSpheretf,
            ValueLocation_CellCentered,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vSpheretf,
            ValueLocation_CellCentered,
            cellVarNames,
            cellVarLocation
        );

        tecplotWriter::getTecplotNames
        (
            vSymmtf,
            ValueLocation_CellCentered,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vSymmtf,
            ValueLocation_CellCentered,
            cellVarNames,
            cellVarLocation
        );

        tecplotWriter::getTecplotNames
        (
            vtf,
            ValueLocation_CellCentered,
            varNames,
            varLocation
        );
        tecplotWriter::getTecplotNames
        (
            vtf,
            ValueLocation_CellCentered,
            cellVarNames,
            cellVarLocation
        );



        // pointFields
        tecplotWriter::getTecplotNames
        (
            psf,
            ValueLocation_Nodal,
            varNames,
            varLocation
        );

        tecplotWriter::getTecplotNames
        (
            pvf,
            ValueLocation_Nodal,
            varNames,
            varLocation
        );

        // strandID (= piece id. Gets incremented for every piece of geometry
        // that is output)
        INTEGER4 strandID = 1;


        if (meshState != polyMesh::UNCHANGED)
        {
            if (doWriteInternal)
            {
                // Output mesh and fields
                fileName vtkFileName
                (
                    fvPath/vtkName
                  + "_"
                  + timeDesc
                  + ".plt"
                );

                tecplotWriter writer(runTime);

                string allVarNames = string("X Y Z ") + varNames;
                DynamicList<INTEGER4> allVarLocation;
                allVarLocation.append(ValueLocation_Nodal);
                allVarLocation.append(ValueLocation_Nodal);
                allVarLocation.append(ValueLocation_Nodal);
                allVarLocation.append(varLocation);

                writer.writeInit
                (
                    runTime.caseName(),
                    allVarNames,
                    vtkFileName,
                    DataFileType_Full
                );

                writer.writePolyhedralZone
                (
                    mesh.name(),        // regionName
                    strandID++,         // strandID
                    mesh,
                    allVarLocation,
                    nFaceNodes
                );

                // Write coordinates
                writer.writeField(mesh.points().component(0)());
                writer.writeField(mesh.points().component(1)());
                writer.writeField(mesh.points().component(2)());

                // Write all fields
                forAll(vsf, i)
                {
                    writer.writeField(vsf[i]);
                }
                forAll(vvf, i)
                {
                    writer.writeField(vvf[i]);
                }
                forAll(vSpheretf, i)
                {
                    writer.writeField(vSpheretf[i]);
                }
                forAll(vSymmtf, i)
                {
                    writer.writeField(vSymmtf[i]);
                }
                forAll(vtf, i)
                {
                    writer.writeField(vtf[i]);
                }

                forAll(psf, i)
                {
                    writer.writeField(psf[i]);
                }
                forAll(pvf, i)
                {
                    writer.writeField(pvf[i]);
                }

                writer.writeConnectivity(mesh);
                writer.writeEnd();
            }
        }
        else
        {
            if (doWriteInternal)
            {
                if (timeI == 0)
                {
                    // Output static mesh only
                    fileName vtkFileName
                    (
                        fvPath/vtkName
                      + "_grid_"
                      + timeDesc
                      + ".plt"
                    );

                    tecplotWriter writer(runTime);

                    writer.writeInit
                    (
                        runTime.caseName(),
                        "X Y Z",
                        vtkFileName,
                        DataFileType_Grid
                    );

                    writer.writePolyhedralZone
                    (
                        mesh.name(),        // regionName
                        strandID,           // strandID
                        mesh,
                        List<INTEGER4>(3, ValueLocation_Nodal),
                        nFaceNodes
                    );

                    // Write coordinates
                    writer.writeField(mesh.points().component(0)());
                    writer.writeField(mesh.points().component(1)());
                    writer.writeField(mesh.points().component(2)());

                    writer.writeConnectivity(mesh);
                    writer.writeEnd();
                }

                // Output solution file
                fileName vtkFileName
                (
                    fvPath/vtkName
                  + "_"
                  + timeDesc
                  + ".plt"
                );

                tecplotWriter writer(runTime);

                writer.writeInit
                (
                    runTime.caseName(),
                    varNames,
                    vtkFileName,
                    DataFileType_Solution
                );

                writer.writePolyhedralZone
                (
                    mesh.name(),        // regionName
                    strandID++,         // strandID
                    mesh,
                    varLocation,
                    0
                );

                // Write all fields
                forAll(vsf, i)
                {
                    writer.writeField(vsf[i]);
                }
                forAll(vvf, i)
                {
                    writer.writeField(vvf[i]);
                }
                forAll(vSpheretf, i)
                {
                    writer.writeField(vSpheretf[i]);
                }
                forAll(vSymmtf, i)
                {
                    writer.writeField(vSymmtf[i]);
                }
                forAll(vtf, i)
                {
                    writer.writeField(vtf[i]);
                }

                forAll(psf, i)
                {
                    writer.writeField(psf[i]);
                }
                forAll(pvf, i)
                {
                    writer.writeField(pvf[i]);
                }
                writer.writeEnd();
            }
        }


        //---------------------------------------------------------------------
        //
        // Write faceSet
        //
        //---------------------------------------------------------------------

        if (args.optionFound("faceSet"))
        {
            // Load the faceSet
            const word setName = args["faceSet"];
            labelList faceLabels(faceSet(mesh, setName).toc());

            // Filename as if patch with same name.
            mkDir(fvPath/setName);

            fileName patchFileName
            (
                fvPath/setName/setName
              + "_"
              + timeDesc
              + ".plt"
            );

            Info<< "    FaceSet   : " << patchFileName << endl;

            tecplotWriter writer(runTime);

            string allVarNames = string("X Y Z ") + cellVarNames;
            DynamicList<INTEGER4> allVarLocation;
            allVarLocation.append(ValueLocation_Nodal);
            allVarLocation.append(ValueLocation_Nodal);
            allVarLocation.append(ValueLocation_Nodal);
            allVarLocation.append(cellVarLocation);

            writer.writeInit
            (
                runTime.caseName(),
                cellVarNames,
                patchFileName,
                DataFileType_Full
            );

            const indirectPrimitivePatch ipp
            (
                IndirectList<face>(mesh.faces(), faceLabels),
                mesh.points()
            );

            writer.writePolygonalZone
            (
                setName,
                strandID++,
                ipp,
                allVarLocation
            );

            // Write coordinates
            writer.writeField(ipp.localPoints().component(0)());
            writer.writeField(ipp.localPoints().component(1)());
            writer.writeField(ipp.localPoints().component(2)());

            // Write all volfields
            forAll(vsf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vsf[i])(),
                        faceLabels
                    )()
                );
            }
            forAll(vvf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vvf[i])(),
                        faceLabels
                    )()
                );
            }
            forAll(vSpheretf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vSpheretf[i])(),
                        faceLabels
                    )()
                );
            }
            forAll(vSymmtf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vSymmtf[i])(),
                        faceLabels
                    )()
                );
            }
            forAll(vtf, i)
            {
                writer.writeField
                (
                    writer.getFaceField
                    (
                        linearInterpolate(vtf[i])(),
                        faceLabels
                    )()
                );
            }
            writer.writeConnectivity(ipp);

            continue;
        }



        //---------------------------------------------------------------------
        //
        // Write patches as multi-zone file
        //
        //---------------------------------------------------------------------

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        labelList patchIDs(getSelectedPatches(patches, excludePatches));

        mkDir(fvPath/"boundaryMesh");

        fileName patchFileName;

        if (vMesh.useSubMesh())
        {
            patchFileName =
                fvPath/"boundaryMesh"/cellSetName
              + "_"
              + timeDesc
              + ".plt";
        }
        else
        {
            patchFileName =
                fvPath/"boundaryMesh"/"boundaryMesh"
              + "_"
              + timeDesc
              + ".plt";
        }

        Info<< "    Combined patches     : " << patchFileName << endl;

        tecplotWriter writer(runTime);

        string allVarNames = string("X Y Z ") + varNames;
        DynamicList<INTEGER4> allVarLocation;
        allVarLocation.append(ValueLocation_Nodal);
        allVarLocation.append(ValueLocation_Nodal);
        allVarLocation.append(ValueLocation_Nodal);
        allVarLocation.append(varLocation);

        writer.writeInit
        (
            runTime.caseName(),
            allVarNames,
            patchFileName,
            DataFileType_Full
        );

        forAll(patchIDs, i)
        {
            label patchID = patchIDs[i];
            const polyPatch& pp = patches[patchID];
            // INTEGER4 strandID = 1 + i;

            if (pp.size() > 0)
            {
                Info<< "    Writing patch " << patchID << "\t" << pp.name()
                    << "\tstrand:" << strandID << nl << endl;

                const indirectPrimitivePatch ipp
                (
                    IndirectList<face>(pp, identity(pp.size())),
                    pp.points()
                );

                writer.writePolygonalZone
                (
                    pp.name(),
                    strandID++,     // strandID,
                    ipp,
                    allVarLocation
                );

                // Write coordinates
                writer.writeField(ipp.localPoints().component(0)());
                writer.writeField(ipp.localPoints().component(1)());
                writer.writeField(ipp.localPoints().component(2)());

                // Write all fields
                forAll(vsf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vsf[i],
                            patchID
                        )()
                    );
                }
                forAll(vvf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vvf[i],
                            patchID
                        )()
                    );
                }
                forAll(vSpheretf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vSpheretf[i],
                            patchID
                        )()
                    );
                }
                forAll(vSymmtf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vSymmtf[i],
                            patchID
                        )()
                    );
                }
                forAll(vtf, i)
                {
                    writer.writeField
                    (
                        writer.getPatchField
                        (
                            nearCellValue,
                            vtf[i],
                            patchID
                        )()
                    );
                }

                forAll(psf, i)
                {
                    writer.writeField
                    (
                        psf[i].boundaryField()[patchID].patchInternalField()()
                    );
                }
                forAll(pvf, i)
                {
                    writer.writeField
                    (
                        pvf[i].boundaryField()[patchID].patchInternalField()()
                    );
                }

                writer.writeConnectivity(ipp);
            }
            else
            {
                Info<< "    Skipping zero sized patch " << patchID
                    << "\t" << pp.name()
                    << nl << endl;
            }
        }
        writer.writeEnd();
        Info<< endl;



        //---------------------------------------------------------------------
        //
        // Write faceZones as multi-zone file
        //
        //---------------------------------------------------------------------

        const faceZoneMesh& zones = mesh.faceZones();

        if (doFaceZones && zones.size() > 0)
        {
            mkDir(fvPath/"faceZoneMesh");

            fileName patchFileName;

            if (vMesh.useSubMesh())
            {
                patchFileName =
                    fvPath/"faceZoneMesh"/cellSetName
                  + "_"
                  + timeDesc
                  + ".plt";
            }
            else
            {
                patchFileName =
                    fvPath/"faceZoneMesh"/"faceZoneMesh"
                  + "_"
                  + timeDesc
                  + ".plt";
            }

            Info<< "    FaceZone  : " << patchFileName << endl;

            tecplotWriter writer(runTime);

            string allVarNames = string("X Y Z ") + cellVarNames;
            DynamicList<INTEGER4> allVarLocation;
            allVarLocation.append(ValueLocation_Nodal);
            allVarLocation.append(ValueLocation_Nodal);
            allVarLocation.append(ValueLocation_Nodal);
            allVarLocation.append(cellVarLocation);

            writer.writeInit
            (
                runTime.caseName(),
                allVarNames,
                patchFileName,
                DataFileType_Full
            );

            forAll(zones, zoneI)
            {
                const faceZone& pp = zones[zoneI];

                if (pp.size() > 0)
                {
                    const indirectPrimitivePatch ipp
                    (
                        IndirectList<face>(mesh.faces(), pp),
                        mesh.points()
                    );

                    writer.writePolygonalZone
                    (
                        pp.name(),
                        strandID++, //1+patchIDs.size()+zoneI,    // strandID,
                        ipp,
                        allVarLocation
                    );

                    // Write coordinates
                    writer.writeField(ipp.localPoints().component(0)());
                    writer.writeField(ipp.localPoints().component(1)());
                    writer.writeField(ipp.localPoints().component(2)());

                    // Write all volfields
                    forAll(vsf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vsf[i])(),
                                pp
                            )()
                        );
                    }
                    forAll(vvf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vvf[i])(),
                                pp
                            )()
                        );
                    }
                    forAll(vSpheretf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vSpheretf[i])(),
                                pp
                            )()
                        );
                    }
                    forAll(vSymmtf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vSymmtf[i])(),
                                pp
                            )()
                        );
                    }
                    forAll(vtf, i)
                    {
                        writer.writeField
                        (
                            writer.getFaceField
                            (
                                linearInterpolate(vtf[i])(),
                                pp
                            )()
                        );
                    }

                    writer.writeConnectivity(ipp);
                }
                else
                {
                    Info<< "    Skipping zero sized faceZone " << zoneI
                        << "\t" << pp.name()
                        << nl << endl;
                }
            }
            writer.writeEnd();
            Info<< endl;
        }



        //---------------------------------------------------------------------
        //
        // Write lagrangian data
        //
        //---------------------------------------------------------------------

        fileNameList cloudDirs
        (
            readDir
            (
                runTime.timePath()/regionPrefix/cloud::prefix,
                fileType::directory
            )
        );

        forAll(cloudDirs, cloudI)
        {
            IOobjectList sprayObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudDirs[cloudI]
            );

            IOobject* positionsPtr = sprayObjs.lookup("positions");

            if (positionsPtr)
            {
                mkDir(fvPath/cloud::prefix/cloudDirs[cloudI]);

                fileName lagrFileName
                (
                    fvPath/cloud::prefix/cloudDirs[cloudI]/cloudDirs[cloudI]
                  + "_" + timeDesc + ".plt"
                );

                Info<< "    Lagrangian: " << lagrFileName << endl;

                wordList labelNames(sprayObjs.names(labelIOField::typeName));
                Info<< "        labels            :";
                print(Info, labelNames);

                wordList scalarNames(sprayObjs.names(scalarIOField::typeName));
                Info<< "        scalars           :";
                print(Info, scalarNames);

                wordList vectorNames(sprayObjs.names(vectorIOField::typeName));
                Info<< "        vectors           :";
                print(Info, vectorNames);

                // wordList sphereNames
                //(
                //    sprayObjs.names
                //    (
                //        sphericalTensorIOField::typeName
                //    )
                //);
                // Info<< "        spherical tensors :";
                // print(Info, sphereNames);
                //
                // wordList symmNames
                //(
                //    sprayObjs.names
                //    (
                //        symmTensorIOField::typeName
                //    )
                //);
                // Info<< "        symm tensors      :";
                // print(Info, symmNames);
                //
                // wordList tensorNames
                //(
                //    sprayObjs.names(tensorIOField::typeName)
                //);
                // Info<< "        tensors           :";
                // print(Info, tensorNames);


                // Load cloud positions
                passiveParticleCloud parcels(mesh, cloudDirs[cloudI]);

                // Get positions as pointField
                pointField positions(parcels.size());
                label n = 0;
                forAllConstIter(passiveParticleCloud, parcels, elmnt)
                {
                    positions[n++] = elmnt().position();
                }


                string allVarNames = string("X Y Z");
                DynamicList<INTEGER4> allVarLocation;
                allVarLocation.append(ValueLocation_Nodal);
                allVarLocation.append(ValueLocation_Nodal);
                allVarLocation.append(ValueLocation_Nodal);

                tecplotWriter::getTecplotNames<label>
                (
                    labelNames,
                    ValueLocation_Nodal,
                    allVarNames,
                    allVarLocation
                );

                tecplotWriter::getTecplotNames<scalar>
                (
                    scalarNames,
                    ValueLocation_Nodal,
                    allVarNames,
                    allVarLocation
                );

                tecplotWriter::getTecplotNames<vector>
                (
                    vectorNames,
                    ValueLocation_Nodal,
                    allVarNames,
                    allVarLocation
                );


                tecplotWriter writer(runTime);

                writer.writeInit
                (
                    runTime.caseName(),
                    allVarNames,
                    lagrFileName,
                    DataFileType_Full
                );

                writer.writeOrderedZone
                (
                    cloudDirs[cloudI],
                    strandID++,     // strandID,
                    parcels.size(),
                    allVarLocation
                );

                // Write coordinates
                writer.writeField(positions.component(0)());
                writer.writeField(positions.component(1)());
                writer.writeField(positions.component(2)());

                // labelFields
                forAll(labelNames, i)
                {
                    IOField<label> fld
                    (
                        IOobject
                        (
                            labelNames[i],
                            mesh.time().timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        )
                    );

                    scalarField sfld(fld.size());
                    forAll(fld, j)
                    {
                        sfld[j] = scalar(fld[j]);
                    }
                    writer.writeField(sfld);
                }
                // scalarFields
                forAll(scalarNames, i)
                {
                    IOField<scalar> fld
                    (
                        IOobject
                        (
                            scalarNames[i],
                            mesh.time().timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        )
                    );
                    writer.writeField(fld);
                }
                // vectorFields
                forAll(vectorNames, i)
                {
                    IOField<vector> fld
                    (
                        IOobject
                        (
                            vectorNames[i],
                            mesh.time().timeName(),
                            cloud::prefix/cloudDirs[cloudI],
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        )
                    );
                    writer.writeField(fld);
                }

                writer.writeEnd();
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
