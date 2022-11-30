/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "mapClouds.H"
#include "fvMeshToFvMesh.H"
#include "IOobjectList.H"
#include "OSspecific.H"
#include "passiveParticleCloud.H"
#include "patchToPatchTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
List<Type> gatherAndFlatten(const List<Type>& l)
{
    List<List<Type>> procLs(Pstream::nProcs());
    procLs[Pstream::myProcNo()] = l;

    Pstream::gatherList(procLs);
    Pstream::scatterList(procLs);

    return
        HashSet<Type>
        (
            ListListOps::combine<List<Type>>
            (
                procLs,
                accessOp<List<Type>>()
            )
        ).sortedToc();
}


template<class ReadIOField, class WriteIOField>
void mapCloudTypeFields
(
    const fileName& cloudDir,
    const IOobjectList& objects,
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const distributionMap& map
)
{
    const wordList fieldNames =
        gatherAndFlatten(objects.lookupClass(ReadIOField::typeName).names());

    forAll(fieldNames, fieldi)
    {
        const word& fieldName = fieldNames[fieldi];

        Info<< "    mapping lagrangian field " << fieldName << endl;

        const IOobject* fieldIOPtr = objects.lookup(fieldName);

        // Read the field
        ReadIOField field
        (
            fieldIOPtr != nullptr
          ? *fieldIOPtr
          : IOobject
            (
                fieldName,
                srcMesh.time().name(),
                cloud::prefix/cloudDir,
                srcMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // Distribute the field
        map.distribute(field);

        // Write it out into the target case
        if (field.size())
        {
            WriteIOField
            (
                IOobject
                (
                    fieldName,
                    tgtMesh.time().name(),
                    cloud::prefix/cloudDir,
                    tgtMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                move(field)
            ).write();
        }
    }
}


template<class Type>
void mapCloudTypeFieldsAndFieldFields
(
    const fileName& cloudDir,
    const IOobjectList& objects,
    const polyMesh& srcMesh,
    const polyMesh& tgtMesh,
    const distributionMap& map
)
{
    mapCloudTypeFields
    <
        IOField<Type>,
        IOField<Type>
    >
    (
        cloudDir,
        objects,
        srcMesh,
        tgtMesh,
        map
    );

    mapCloudTypeFields
    <
        IOField<Field<Type>>,
        CompactIOField<Field<Type>>
    >
    (
        cloudDir,
        objects,
        srcMesh,
        tgtMesh,
        map
    );

    mapCloudTypeFields
    <
        CompactIOField<Field<Type>>,
        CompactIOField<Field<Type>>
    >
    (
        cloudDir,
        objects,
        srcMesh,
        tgtMesh,
        map
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::mapClouds(const fvMeshToFvMesh& interp)
{
    const polyMesh& srcMesh = interp.srcMesh();
    const polyMesh& tgtMesh = interp.tgtMesh();

    // Determine the clouds present in this mesh
    const fileNameList cloudDirs =
        gatherAndFlatten
        (
            readDir
            (
                srcMesh.time().timePath()/cloud::prefix,
                fileType::directory
            )
        );

    forAll(cloudDirs, cloudi)
    {
        Info<< nl << "Mapping cloud " << cloudDirs[cloudi] << endl;

        // Read the source cloud positions
        passiveParticleCloud srcCloud(srcMesh, cloudDirs[cloudi], false);

        Info<< "    read "
            << returnReduce(srcCloud.size(), sumOp<label>())
            << " parcels from source mesh." << endl;

        // Unpack into position and cell lists and construct a send map
        pointField positions(srcCloud.size());
        labelList tgtCells(srcCloud.size(), -1);
        List<DynamicList<label>> sendsDyn(Pstream::nProcs());
        {
            label srcParticlei = 0;
            forAllConstIter(Cloud<passiveParticle>, srcCloud, iter)
            {
                const point pos = iter().position(srcMesh);

                const remote tgtProcCell =
                    interp.srcToTgtPoint(iter().cell(), pos);

                positions[srcParticlei] = pos;
                tgtCells[srcParticlei] = tgtProcCell.elementi;

                sendsDyn[tgtProcCell.proci].append(srcParticlei ++);
            }
        }
        labelListList sends;
        patchToPatchTools::transferListList(sends, sendsDyn);

        // Build a distribution map from the send map
        autoPtr<distributionMap> mapPtr =
            patchToPatchTools::constructDistributionMap(sends);
        const distributionMap& map = mapPtr();

        // Distribute the positions and target cell indices
        map.distribute(positions);
        map.distribute(tgtCells);

        // Construct the target cloud positions
        passiveParticleCloud tgtCloud
        (
            tgtMesh,
            cloudDirs[cloudi],
            IDLList<passiveParticle>()
        );
        forAll(positions, tgtParticlei)
        {
            tgtCloud.addParticle
            (
                new passiveParticle
                (
                    tgtMesh,
                    positions[tgtParticlei],
                    tgtCells[tgtParticlei]
                )
            );
        }

        Info<< "    mapped "
            << returnReduce(tgtCloud.size(), sumOp<label>())
            << " parcels to the target mesh." << endl;

        // Write the positions
        IOPosition<passiveParticleCloud>(tgtCloud).write();

        // Search for Lagrangian objects for this time
        IOobjectList objects
        (
            srcMesh,
            srcMesh.time().name(),
            cloud::prefix/cloudDirs[cloudi]
        );

        // Map and write the fields
        #define MapCloudTypeFields(Type, nullArg)                          \
            mapCloudTypeFieldsAndFieldFields<Type>                         \
            (                                                              \
                cloudDirs[cloudi],                                         \
                objects,                                                   \
                srcMesh,                                                   \
                tgtMesh,                                                   \
                map                                                        \
            );
        MapCloudTypeFields(label, );
        FOR_ALL_FIELD_TYPES(MapCloudTypeFields);
        #undef MapCloudTypeFields
    }
}


// ************************************************************************* //
