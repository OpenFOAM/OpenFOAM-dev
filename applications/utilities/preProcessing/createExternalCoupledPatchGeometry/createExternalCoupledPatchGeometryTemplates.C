/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "createExternalCoupledPatchGeometryTemplates.H"
#include "externalCoupledMixedFvPatchField.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void processField
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const word& fieldName,
    label& processed
)
{
    if (processed != -1)
    {
        return;
    }

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const word timeName(mesh.time().timeName());

    IOobjectList fieldObjbjects(objects.lookupClass(fieldType::typeName));

    if (fieldObjbjects.lookup(fieldName) != nullptr)
    {
        fieldType vtf(*fieldObjbjects.lookup(fieldName), mesh);
        const typename fieldType::Boundary& bf =
            vtf.boundaryField();

        forAll(bf, patchi)
        {
            if (isA<externalCoupledMixedFvPatchField<Type>>(bf[patchi]))
            {
                Info<< "Generating external coupled geometry for field "
                    << fieldName << endl;

                const externalCoupledMixedFvPatchField<Type>& pf =
                    refCast<const externalCoupledMixedFvPatchField<Type>>
                    (
                        bf[patchi]
                    );

                pf.writeGeometry();
                processed = 1;

                break;
            }
        }

        if (processed != 1)
        {
            processed = 0;

            Info<< "Field " << fieldName << " found, but does not have any "
                << externalCoupledMixedFvPatchField<Type>::typeName
                << " boundary conditions" << endl;
        }
    }
}


// ************************************************************************* //
