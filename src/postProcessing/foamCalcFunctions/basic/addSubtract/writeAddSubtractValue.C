/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
void Foam::calcTypes::addSubtract::writeAddSubtractValue
(
    const IOobject& baseHeader,
    const string& valueStr,
    const fvMesh& mesh,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (baseHeader.headerClassName() == fieldType::typeName)
    {
        if (resultName_ == "")
        {
            if (calcMode_ == ADD)
            {
                resultName_ = baseHeader.name() + "_add_value";
            }
            else
            {
                resultName_ = baseHeader.name() + "_subtract_value";
            }
        }

        Type value;
        IStringStream(valueStr)() >> value;

        Info<< "    Reading " << baseHeader.name() << endl;
        fieldType baseField(baseHeader, mesh);

        fieldType newField
        (
            IOobject
            (
                resultName_,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ
            ),
            baseField
        );

        Info<< "    Calculating " << resultName_ << endl;
        if (calcMode_ == ADD)
        {
            newField == baseField
                + dimensioned<Type>("value", baseField.dimensions(), value);
        }
        else
        {
            newField == baseField
                - dimensioned<Type>("value", baseField.dimensions(), value);
        }

        newField.write();

        processed = true;
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
