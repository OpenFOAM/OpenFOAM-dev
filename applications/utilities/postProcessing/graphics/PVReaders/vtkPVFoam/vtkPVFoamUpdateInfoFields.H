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

InClass
    vtkPVFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkPVFoamUpdateInfoFields_H
#define vtkPVFoamUpdateInfoFields_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class patchType, class meshType>
void Foam::vtkPVFoam::updateInfoFields
(
    vtkDataArraySelection* select
)
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPVFoam::updateInfoFields <"
            << meshType::Mesh::typeName
            << "> [meshPtr=" << (meshPtr_ ? "set" : "nullptr") << "]"
            << endl;
    }

    // Use the db directly since this might be called without a mesh,
    // but the region must get added back in
    word regionPrefix;
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        regionPrefix = meshRegion_;
    }

    const instantList times = dbPtr_().times();

    stringList enabledEntries;

    if (select->GetNumberOfArrays() == 0 && !meshPtr_)
    {
        const wordReList defaultFieldRes
        (
            configDict_.lookupOrDefault
            (
                "defaultFields",
                wordReList(wordList{"p", "U"})
            )
        );

        wordHashSet objectNameSet;
        forAll(times, timei)
        {
            // Search for list of objects for this time and mesh region
            IOobjectList objects(dbPtr_(), times[timei].name(), regionPrefix);

            forAllConstIter(IOobjectList, objects, iter)
            {
                objectNameSet.insert(iter.key());
            }
        }

        const wordList objectNames(objectNameSet.toc());
        const labelList defaultFields
        (
            findStrings(defaultFieldRes, objectNames)
        );

        enabledEntries.setSize(defaultFields.size());
        forAll(defaultFields, i)
        {
            enabledEntries[i] = objectNames[defaultFields[i]];
        }
    }
    else
    {
        // Preserve the enabled selections
        enabledEntries = getSelectedArrayEntries(select);
    }

    select->RemoveAllArrays();

    forAll(times, timei)
    {
        // Search for list of objects for this time and mesh region
        IOobjectList objects(dbPtr_(), times[timei].name(), regionPrefix);

        //- Add volume fields to GUI
        addToSelection<GeometricField<scalar, patchType, meshType>>
        (
            select,
            objects
        );
        addToSelection<GeometricField<vector, patchType, meshType>>
        (
            select,
            objects
        );
        addToSelection<GeometricField<sphericalTensor, patchType, meshType>>
        (
            select,
            objects
        );
        addToSelection<GeometricField<symmTensor, patchType, meshType>>
        (
            select,
            objects
        );
        addToSelection<GeometricField<tensor, patchType, meshType>>
        (
            select,
            objects
        );
    }

    // Restore the enabled selections
    setSelectedArrayEntries(select, enabledEntries);

    if (debug)
    {
        Info<< "<end> Foam::vtkPVFoam::updateInfoFields" << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
