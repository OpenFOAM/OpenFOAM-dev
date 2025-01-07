/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "toSubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class DimensionedFieldToDimensionedSubField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class DimensionedFieldToDimensionedSubField
:
    public DimensionedField<Type, GeoMesh, SubField>
{
    // Private Data

        //- The primitive field
        autoPtr<Field<Type>> primitiveFieldPtr_;


public:

    // Constructors

        //- Construct from a temporary field
        DimensionedFieldToDimensionedSubField
        (
            const tmp<DimensionedField<Type, GeoMesh, Field>>& tField
        )
        :
            DimensionedField<Type, GeoMesh, SubField>
            (
                tField(),
                tField->mesh(),
                tField->dimensions(),
                SubField<Type>::null()
            ),
            primitiveFieldPtr_(nullptr)
        {
            if (tField.isTmp())
            {
                primitiveFieldPtr_ =
                    new Field<Type>(tField.ref().primitiveFieldRef(), true);

                this->shallowCopy(primitiveFieldPtr_());
            }
            else
            {
                this->shallowCopy(tField->primitiveField());
            }

            tField.clear();
        }


    //- Destructor
    virtual ~DimensionedFieldToDimensionedSubField()
    {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, Foam::SubField>>
Foam::toSubField(const DimensionedField<Type, GeoMesh, Field>& f)
{
    return
        tmp<DimensionedField<Type, GeoMesh, SubField>>
        (
            // Use the converter class to deep-copy the IOobject and
            // dimensionSet and shallow-copy the Field
            /*
            new DimensionedFieldToDimensionedSubField<Type, GeoMesh>
            (
                tmp<DimensionedField<Type, GeoMesh, Field>>(f)
            )
            */

            // Cast. The layout of Field and SubField are the same, and because
            // this is a const-reference tmp, we don't have to worry about the
            // SubField not destructing the data. This prevents the converter
            // class' copy of the IOobject and dimensionSet.
            reinterpret_cast<const DimensionedField<Type, GeoMesh, SubField>&>
            (
                f
            )
        );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, Foam::SubField>>
Foam::toSubField(const DimensionedField<Type, GeoMesh, SubField>& f)
{
    return tmp<DimensionedField<Type, GeoMesh, SubField>>(f);
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, Foam::SubField>>
Foam::toSubField(const tmp<DimensionedField<Type, GeoMesh, Field>>& tf)
{
    return
        tmp<DimensionedField<Type, GeoMesh, SubField>>
        (
            new DimensionedFieldToDimensionedSubField<Type, GeoMesh>(tf),
            true
        );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, Foam::SubField>>
Foam::toSubField(const tmp<DimensionedField<Type, GeoMesh, SubField>>& tf)
{
    return tmp<DimensionedField<Type, GeoMesh, SubField>>(tf, true);
}


template<class Type, class GeoMesh, class ... Args>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh, Foam::SubField>>
Foam::toSubField(const Args& ... args)
{
    return toSubField(DimensionedField<Type, GeoMesh, Field>::New(args ...));
}


// ************************************************************************* //
