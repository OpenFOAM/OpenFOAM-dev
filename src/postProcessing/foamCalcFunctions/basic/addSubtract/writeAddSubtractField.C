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

\*---------------------------------------------------------------------------*/

template<class Type>
void Foam::calcTypes::addSubtract::writeAddSubtractField
(
    const IOobject& baseHeader,
    const IOobject& addHeader,
    const fvMesh& mesh,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if
    (
        baseHeader.headerClassName() == fieldType::typeName
     && baseHeader.headerClassName() == addHeader.headerClassName()
    )
    {
        if (resultName_ == "")
        {
            if (calcMode_ == ADD)
            {
                resultName_ = baseHeader.name() + "_add_" + addHeader.name();
            }
            else
            {
                resultName_ = baseHeader.name() + "_subtract_"
                    + addHeader.name();
            }
        }

        Info<< "    Reading " << baseHeader.name() << endl;
        fieldType baseField(baseHeader, mesh);

        Info<< "    Reading " << addHeader.name() << endl;
        fieldType addField(addHeader, mesh);

        if (baseField.dimensions() == addField.dimensions())
        {
            Info<< "    Calculating " << resultName_ << endl;

            fieldType newField
            (
                IOobject
                (
                    resultName_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                calcMode_ == ADD
              ? (baseField + addField)()
              : (baseField - addField)()
            );
            newField.write();
        }
        else
        {
            Info<< "    Cannot calculate " << resultName_ << nl
                << "    - inconsistent dimensions: "
                << baseField.dimensions() << " - " << addField.dimensions()
                << endl;
        }

        processed = true;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
