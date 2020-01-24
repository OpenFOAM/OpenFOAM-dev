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

Application
    foamToTetDualMesh

Description
    Converts polyMesh results to tetDualMesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "Time.H"
#include "IOobjectList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ReadGeoField, class MappedGeoField>
void ReadAndMapFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const fvMesh& tetDualMesh,
    const labelList& map,
    const typename MappedGeoField::value_type& nullValue,
    PtrList<MappedGeoField>& tetFields
)
{
    typedef typename MappedGeoField::value_type Type;

    // Search list of objects for wanted type
    IOobjectList fieldObjects(objects.lookupClass(ReadGeoField::typeName));

    tetFields.setSize(fieldObjects.size());

    label i = 0;
    forAllConstIter(IOobjectList, fieldObjects, iter)
    {
        Info<< "Converting " << ReadGeoField::typeName << ' ' << iter.key()
            << endl;

        ReadGeoField readField(*iter(), mesh);

        tetFields.set
        (
            i,
            new MappedGeoField
            (
                IOobject
                (
                    readField.name(),
                    readField.instance(),
                    readField.local(),
                    tetDualMesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    readField.registerObject()
                ),
                pointMesh::New(tetDualMesh),
                dimensioned<Type>
                (
                    "zero",
                    readField.dimensions(),
                    Zero
                )
            )
        );

        Field<Type>& fld = tetFields[i].primitiveFieldRef();

        // Map from read field. Set unmapped entries to nullValue.
        fld.setSize(map.size(), nullValue);
        forAll(map, pointi)
        {
            label index = map[pointi];

            if (index > 0)
            {
                label celli = index-1;
                fld[pointi] = readField[celli];
            }
            else if (index < 0)
            {
                label facei = -index-1;
                label bFacei = facei - mesh.nInternalFaces();
                if (bFacei >= 0)
                {
                    label patchi = mesh.boundaryMesh().patchID()[bFacei];
                    label localFacei = mesh.boundaryMesh()[patchi].whichFace
                    (
                        facei
                    );
                    fld[pointi] = readField.boundaryField()[patchi][localFacei];
                }
                // else
                //{
                //    FatalErrorInFunction
                //        << "Face " << facei << " from index " << index
                //        << " is not a boundary face." << abort(FatalError);
                //}

            }
            // else
            //{
            //    WarningInFunction
            //        << "Point " << pointi << " at "
            //        << tetDualMesh.points()[pointi]
            //        << " has no dual correspondence." << endl;
            //}
        }
        tetFields[i].correctBoundaryConditions();

        i++;
    }
}




int main(int argc, char *argv[])
{
    #include "addOverwriteOption.H"
    #include "addTimeOptions.H"

    #include "setRootCase.H"
    #include "createTime.H"
    // Get times list
    instantList Times = runTime.times();
    #include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);


    // Read the mesh
    #include "createMesh.H"

    // Read the tetDualMesh
    Info<< "Create tetDualMesh for time = "
        << runTime.timeName() << nl << endl;

    fvMesh tetDualMesh
    (
        IOobject
        (
            "tetDualMesh",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
    // From tet vertices to poly cells/faces
    const labelIOList pointDualAddressing
    (
        IOobject
        (
            "pointDualAddressing",
            tetDualMesh.facesInstance(),
            tetDualMesh.meshSubDir,
            tetDualMesh,
            IOobject::MUST_READ
        )
    );

    if (pointDualAddressing.size() != tetDualMesh.nPoints())
    {
            FatalErrorInFunction
                << "Size " << pointDualAddressing.size()
                << " of addressing map "
                << pointDualAddressing.localObjectPath()
                << " differs from number of points in mesh "
                << tetDualMesh.nPoints()
                << exit(FatalError);
    }


    // Some stats on addressing
    label nCells = 0;
    label nPatchFaces = 0;
    label nUnmapped = 0;
    forAll(pointDualAddressing, pointi)
    {
        label index = pointDualAddressing[pointi];

        if (index > 0)
        {
            nCells++;
        }
        else if (index == 0)
        {
            nUnmapped++;
        }
        else
        {
            label facei = -index-1;
            if (facei < mesh.nInternalFaces())
            {
                FatalErrorInFunction
                    << "Face " << facei << " from index " << index
                    << " is not a boundary face."
                    << " nInternalFaces:" << mesh.nInternalFaces()
                    << exit(FatalError);
            }
            else
            {
                nPatchFaces++;
            }
        }
    }

    reduce(nCells, sumOp<label>());
    reduce(nPatchFaces, sumOp<label>());
    reduce(nUnmapped, sumOp<label>());
    Info<< "tetDualMesh points : " << tetDualMesh.nPoints()
        << " of which mapped to" << nl
        << "    cells       : " << nCells << nl
        << "    patch faces : " << nPatchFaces << nl
        << "    not mapped  : " << nUnmapped << nl
        << endl;


    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields, interpolate onto tet points
    PtrList<pointScalarField> psFlds;
    ReadAndMapFields<volScalarField, pointScalarField>
    (
        mesh,
        objects,
        tetDualMesh,
        pointDualAddressing,
        Zero,  // nullValue
        psFlds
    );

    PtrList<pointVectorField> pvFlds;
    ReadAndMapFields<volVectorField, pointVectorField>
    (
        mesh,
        objects,
        tetDualMesh,
        pointDualAddressing,
        Zero,  // nullValue
        pvFlds
    );

    PtrList<pointSphericalTensorField> pstFlds;
    ReadAndMapFields<volSphericalTensorField, pointSphericalTensorField>
    (
        mesh,
        objects,
        tetDualMesh,
        pointDualAddressing,
        Zero,  // nullValue
        pstFlds
    );

    PtrList<pointSymmTensorField> psymmtFlds;
    ReadAndMapFields<volSymmTensorField, pointSymmTensorField>
    (
        mesh,
        objects,
        tetDualMesh,
        pointDualAddressing,
        Zero,  // nullValue
        psymmtFlds
    );

    PtrList<pointTensorField> ptFlds;
    ReadAndMapFields<volTensorField, pointTensorField>
    (
        mesh,
        objects,
        tetDualMesh,
        pointDualAddressing,
        Zero,  // nullValue
        ptFlds
    );

    tetDualMesh.objectRegistry::write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
