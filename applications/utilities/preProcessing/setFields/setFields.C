/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Initialises fields with default and zone values

    The default and zone values are read from a dictionary which defaults to
    system/setFieldsDict, cellZones are used to specify the internal field
    values and faceZone patch field values or by extrapolation from the
    internal field.

Usage
    Any number of fields can be initialised on any number of zones, for
    example the 1D shock tube is initialised by
    \verbatim
        defaultValues
        {
            U       (0 0 0);
            T       348.432;
            p       100000;
        }

        zones
        {
            lowPressure
            {
                type        box;
                zoneType    cell;

                box         (0 -1 -1) (5 1 1);

                values
                {
                    T       278.746;
                    p       10000;
                }
            }
        }
    \endverbatim
    and the water in the tank of the rotatingCube VoF case is initialised by
    \verbatim
        defaultValues
        {
            alpha.water 0;
        }

        zones
        {
            cells
            {
                type        box;
                zoneType    cell;

                box (-1e300 -1e300 -1e300) (1e300 0 1e300);

                values
                {
                    alpha.water 1;
                }
            }

            extrapolatePatches
            {
                "inlet|outlet"   (alpha.water);
            }
        }
    \endverbatim
    which sets the internal values of phase-fraction field and inlet and outlet
    patch values by extrapolation from the internal to provide consistent
    boundary distribution.

    Alternatively the inlet and outlet patch values can be set explicitly on a
    faceZone constructed from the corresponding patches:
    \verbatim
        defaultValues
        {
            alpha.water 0;
        }

        zones
        {
            cells
            {
                type        box;
                zoneType    cell;

                box (-1e300 -1e300 -1e300) (1e300 0 1e300);

                values
                {
                    alpha.water 1;
                }
            }

            patchFaces
            {
                type        box;
                zoneType    face;

                box (-1e300 -1e300 -1e300) (1e300 0 1e300);

                zone
                {
                    type        patch;
                    patches     (inlet outlet);
                }

                values
                {
                    alpha.water 1;
                }
            }
        }
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "zoneGenerator.H"
#include "systemDict.H"

// Backward compatibility
#include "topoSetSource.H"
#include "cellSet.H"
#include "faceSet.H"
#include "setCellField.H"

using namespace Foam;

void setVolFields
(
    const fvMesh& mesh,
    const dictionary& fieldsDict,
    const labelList& selectedCells,
    const PtrDictionary<wordReList>& extrapolatePatches
);

void setPatchFields
(
    const fvMesh& mesh,
    const dictionary& fieldsDict,
    const labelList& selectedCells
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    timeSelector::select0(runTime, args);
    #include "createRegionMeshNoChangers.H"

    const dictionary setFieldsDict(systemDict("setFieldsDict", args, mesh));

    PtrDictionary<wordReList> extrapolatePatches;

    if (setFieldsDict.isDict("extrapolatePatches"))
    {
        const dictionary& extrapolatePatchesDict =
            setFieldsDict.subDict("extrapolatePatches");

        forAll(mesh.boundary(), patchi)
        {
            if (extrapolatePatchesDict.found(mesh.boundary()[patchi].name()))
            {
                extrapolatePatches.insert
                (
                    mesh.boundary()[patchi].name(),
                    new wordReList
                    (
                        extrapolatePatchesDict.lookup
                        (
                            mesh.boundary()[patchi].name()
                        )
                    )
                );
            }
        }
    }

    if (setFieldsDict.found("defaultValues"))
    {
        Info<< "Setting field default values" << endl;
        setVolFields
        (
            mesh,
            setFieldsDict.subDict("defaultValues"),
            labelList::null(),
            extrapolatePatches
        );
        Info<< endl;
    }
    else if (setFieldsDict.found("defaultFieldValues"))
    {
        Info<< "Setting field default values" << nl
            << "    The 'defaultFieldValues' entry is deprecated, "
               "please use 'default'" << endl;

        PtrList<setCellField> defaultFieldValues
        (
            setFieldsDict.lookup("defaultFieldValues"),
            setCellField::iNew(mesh, labelList::null())
        );
        Info<< endl;
    }

    if (setFieldsDict.found("zones"))
    {
        Info<< "Setting field zone values" << endl;

        const dictionary& zonesDict = setFieldsDict.subDict("zones");

        forAllConstIter(dictionary, zonesDict, iter)
        {
            Info<< "Zone: " << iter().keyword() << endl;

            const dictionary& zoneDict = iter().dict();

            autoPtr<zoneGenerator> zg;

            if (zoneDict.found("zoneType"))
            {
                zg = zoneGenerator::New(iter().keyword(), mesh, zoneDict);
            }
            else
            {
                zg = zoneGenerator::New
                (
                    iter().keyword(),
                    zoneTypes::cell,
                    mesh,
                    zoneDict
                );
            }

            const zoneSet zs(zg->generate());

            if (zs.cValid())
            {
                setVolFields
                (
                    mesh,
                    zoneDict.subDict("values"),
                    zs.cZone(),
                    extrapolatePatches
                );
            }

            if (zs.fValid())
            {
                setPatchFields(mesh, zoneDict.subDict("values"), zs.fZone());
            }
        }
    }

    if (setFieldsDict.found("regions"))
    {
        PtrList<entry> regions(setFieldsDict.lookup("regions"));

        Info<< "Setting field region values" << nl
            << "    The 'regions' entry is deprecated, "
               "please use 'zones'" << endl;

        forAll(regions, ri)
        {
            const entry& region = regions[ri];

            autoPtr<topoSetSource> source =
                topoSetSource::New(region.keyword(), mesh, region.dict());

            if (source().setType() == topoSetSource::CELLSETSOURCE)
            {
                cellSet selectedCellSet
                (
                    mesh,
                    "cellSet",
                    mesh.nCells()/10+1  // Reasonable size estimate.
                );

                source->applyToSet
                (
                    topoSetSource::NEW,
                    selectedCellSet
                );

                PtrList<setCellField> fieldValues
                (
                    region.dict().lookup("fieldValues"),
                    setCellField::iNew(mesh, selectedCellSet.toc())
                );
            }
            else if (source().setType() == topoSetSource::FACESETSOURCE)
            {
                faceSet selectedFaceSet
                (
                    mesh,
                    "faceSet",
                    (mesh.nFaces()-mesh.nInternalFaces())/10+1
                );

                source->applyToSet
                (
                    topoSetSource::NEW,
                    selectedFaceSet
                );

                PtrList<setFaceField> fieldValues
                (
                    region.dict().lookup("fieldValues"),
                    setFaceField::iNew(mesh, selectedFaceSet.toc())
                );
            }
        }
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
