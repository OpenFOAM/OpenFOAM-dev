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

#ifndef MapVolFields_H
#define MapVolFields_H

#include "GeometricField.H"
#include "meshToMesh0.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void MapVolFields
(
    const IOobjectList& objects,
    const meshToMesh0& meshToMesh0Interp,
    const meshToMesh0::order& mapOrder
)
{
    const fvMesh& meshSource = meshToMesh0Interp.fromMesh();
    const fvMesh& meshTarget = meshToMesh0Interp.toMesh();

    IOobjectList fields = objects.lookupClass(VolField<Type>::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        typeIOobject<VolField<Type>> fieldTargetIOobject
        (
            fieldIter()->name(),
            meshTarget.time().name(),
            meshTarget,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        if (fieldTargetIOobject.headerOk())
        {
            Info<< "    interpolating " << fieldIter()->name() << endl;

            // Read field fieldSource
            VolField<Type> fieldSource(*fieldIter(), meshSource);

            // Read fieldTarget
            VolField<Type> fieldTarget
            (
                fieldTargetIOobject,
                meshTarget
            );

            // Interpolate field
            meshToMesh0Interp.interpolate
            (
                fieldTarget,
                fieldSource,
                mapOrder
            );

            // Write field
            fieldTarget.write();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
