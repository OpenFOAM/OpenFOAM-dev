/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    dsmcFieldsCalc

Description
    Calculate intensive fields (U and T) from averaged extensive fields from a
    DSMC calculation.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "dsmcCloud.H"
#include "dsmcFields.H"
#include "IOobjectList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<class Type>
    bool addFieldsToList
    (
        const fvMesh& mesh,
        PtrList<GeometricField<Type, fvPatchField, volMesh> >& list,
        const wordList& fieldNames
    )
    {
        typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

        label index = 0;
        forAll(fieldNames, i)
        {
            IOobject obj
            (
                fieldNames[i],
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
            );

            if (obj.headerOk() && obj.headerClassName() == fieldType::typeName)
            {
                list.set(index++, new fieldType(obj, mesh));
            }
            else
            {
                Info<< "Could not find " << fieldNames[i] << endl;

                return false;
            }
        }

        return true;
    }
}


void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    wordList extensiveVSFNames
    (
        IStringStream
        (
            "( \
                rhoNMean \
                rhoMMean \
                linearKEMean \
                internalEMean \
                iDofMean \
            )"
        )()
    );

    PtrList<volScalarField> extensiveVSFs(extensiveVSFNames.size());

    if
    (
        !addFieldsToList
        (
            mesh,
            extensiveVSFs,
            extensiveVSFNames
        )
    )
    {
        return;
    }

    wordList extensiveVVFNames
    (
        IStringStream
        (
            "( \
                momentumMean \
                fDMean \
            )"
        )()
    );

    PtrList<volVectorField> extensiveVVFs(extensiveVVFNames.size());

    if
    (
        !addFieldsToList
        (
            mesh,
            extensiveVVFs,
            extensiveVVFNames
        )
    )
    {
        return;
    }

    dsmcFields dF
    (
        "dsmcFieldsUtility",
        mesh,
        dictionary::null,
        false
    );

    if (writeResults)
    {
        dF.write();
    }
}

// ************************************************************************* //
