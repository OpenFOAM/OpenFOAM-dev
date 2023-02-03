/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "mapGeometricFields.H"
#include "fvMeshToFvMesh.H"
#include "surfaceMesh.H"
#include "pointMesh.H"
#include "IOobjectList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void mapVolTypeFields
(
    const fvMeshToFvMesh& interp,
    const wordReList& cuttingPatches,
    const HashSet<word>& selectedFields,
    const IOobjectList& objects
)
{
    const fvMesh& srcMesh = static_cast<const fvMesh&>(interp.srcMesh());
    const fvMesh& tgtMesh = static_cast<const fvMesh&>(interp.tgtMesh());

    IOobjectList fields = objects.lookupClass(VolField<Type>::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        const word& fieldName = fieldIter()->name();

        if (!selectedFields.empty() && !selectedFields.found(fieldName))
        {
            continue;
        }

        const VolField<Type> fieldSource(*fieldIter(), srcMesh);

        typeIOobject<VolField<Type>> targetIO
        (
            fieldName,
            tgtMesh.time().name(),
            tgtMesh,
            IOobject::READ_IF_PRESENT
        );

        // Warnings about inconsistent execution
        if (targetIO.headerOk() && interp.consistent())
        {
            WarningInFunction
                << "Mapping of field " << fieldName << " will not utilise "
                << "the corresponding field in the target case, as the map is "
                << "consistent (i.e., all patches are mapped)" << endl;
        }
        if (!targetIO.headerOk() && !interp.consistent())
        {
            WarningInFunction
                << "Cannot map field " << fieldName << " because the "
                << "map is not consistent (i.e., not all patches are "
                << "mapped), and there is no corresponding field in "
                << "the target case" << endl;
            continue;
        }
        if (!targetIO.headerOk() && !cuttingPatches.empty())
        {
            WarningInFunction
                << "Cutting patches will not be used for field " << fieldName
                << " because no there is no corresponding field in the target "
                << "case" << endl;
        }

        if (targetIO.headerOk())
        {
            Info<< "    mapping into existing field " << fieldName << endl;

            VolField<Type> fieldTarget(targetIO, tgtMesh);

            fieldTarget.reset
            (
                interp.srcToTgt(fieldSource, fieldTarget, cuttingPatches)
            );

            fieldTarget.write();
        }
        else
        {
            Info<< "    creating new field " << fieldName << endl;

            VolField<Type>(targetIO, interp.srcToTgt(fieldSource)).write();
        }
    }
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::mapGeometricFields
(
    const fvMeshToFvMesh& interp,
    const wordReList& cuttingPatches,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    // Search for list of source objects for this time
    const polyMesh& srcMesh = interp.srcMesh();
    IOobjectList objects(srcMesh, srcMesh.time().name());

    // Map the fields
    #define MapVolTypeFields(Type, nullArg)                                    \
        mapVolTypeFields<Type>                                                 \
        (                                                                      \
            interp,                                                            \
            cuttingPatches,                                                    \
            selectedFields,                                                    \
            objects                                                            \
        );
    FOR_ALL_FIELD_TYPES(MapVolTypeFields);
    #undef MapVolTypeFields
}


// ************************************************************************* //
