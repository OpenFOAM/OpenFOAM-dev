/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "Zonal_DimensionedFvPatchFieldFunction.H"
#include "DimensionedField.H"
#include "zoneGenerator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::Zonal<DimensionedFieldType>::
Zonal
(
    const dictionary& dict,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dict, field),
    value_(dict.lookup<Type>("defaultValue", this->field_.dimensions())),
    zonesDict_(dict.subDict("zones"))
{
    evaluate();
}


template<class DimensionedFieldType>
Foam::DimensionedFieldFunctions::Zonal<DimensionedFieldType>::
Zonal
(
    const Zonal& dff,
    DimensionedFieldType& field
)
:
    DimensionedFieldFunction<DimensionedFieldType>(dff, field),
    value_(dff.value_),
    zonesDict_(dff.zonesDict_)
{}


template<class DimensionedFieldType>
Foam::autoPtr<Foam::DimensionedFieldFunction<DimensionedFieldType>>
Foam::DimensionedFieldFunctions::Zonal<DimensionedFieldType>::clone
(
    DimensionedFieldType& field
) const
{
    return autoPtr<DimensionedFieldFunction<DimensionedFieldType>>
    (
        new Zonal<DimensionedFieldType>(*this, field)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Zonal<DimensionedFieldType>::evaluate()
{
    DimensionedFieldType& field = this->field_;

    const fvMesh& mesh(field.mesh()());

    field.primitiveFieldRef() = value_;

    forAllConstIter(dictionary, zonesDict_, iter)
    {
        const dictionary& zoneDict = iter().dict();

        autoPtr<zoneGenerator> zg
        (
            zoneGenerator::New
            (
                iter().keyword(),
                zoneType<Zone>(),
                mesh,
                zoneDict
            )
        );

        const zoneSet zs(zg->generate());

        if (zs.valid<Zone>())
        {
            const labelList& selected = zs.zone<Zone>();

            const Type value
            (
                zoneDict.lookup<Type>("value", field.dimensions())
            );

            if (&selected == &labelList::null())
            {
                field.primitiveFieldRef() = value;
            }
            else
            {
                forAll(selected, i)
                {
                    const label patchi =
                        mesh.boundaryMesh().whichPatch(selected[i]);
                    if (patchi == field.mesh().index())
                    {
                        field[selected[i] - field.mesh().start()] = value;
                    }
                }
            }
        }
    }
}


template<class DimensionedFieldType>
void Foam::DimensionedFieldFunctions::Zonal<DimensionedFieldType>::write
(
    Ostream& os
) const
{
    writeEntry(os, "defaultValue", this->field_.dimensions(), value_);
    writeEntry(os, "zones", zonesDict_);
}


// ************************************************************************* //
