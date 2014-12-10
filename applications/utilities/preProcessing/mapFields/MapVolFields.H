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

\*---------------------------------------------------------------------------*/

#ifndef MapConsistentVolFields_H
#define MapConsistentVolFields_H

#include "GeometricField.H"
#include "meshToMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class CombineOp>
void MapVolFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    const meshToMesh& interp,
    const CombineOp& cop
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const fvMesh& meshSource = static_cast<const fvMesh&>(interp.srcRegion());
    const fvMesh& meshTarget = static_cast<const fvMesh&>(interp.tgtRegion());

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        const word& fieldName = fieldIter()->name();

        if (selectedFields.empty() || selectedFields.found(fieldName))
        {
            Info<< "    interpolating " << fieldName << endl;

            const fieldType fieldSource(*fieldIter(), meshSource);

            IOobject targetIO
            (
                fieldName,
                meshTarget.time().timeName(),
                meshTarget,
                IOobject::MUST_READ
            );

            if (targetIO.headerOk())
            {
                fieldType fieldTarget(targetIO, meshTarget);

                interp.mapSrcToTgt(fieldSource, cop, fieldTarget);

                fieldTarget.write();
            }
            else
            {
                targetIO.readOpt() = IOobject::NO_READ;

                tmp<fieldType>
                    tfieldTarget(interp.mapSrcToTgt(fieldSource, cop));

                fieldType fieldTarget(targetIO, tfieldTarget);

                fieldTarget.write();
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
