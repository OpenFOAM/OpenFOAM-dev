/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#ifndef DimensionedFieldReuseFunctions_H
#define DimensionedFieldReuseFunctions_H

#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TypeR, class GeoMesh>
tmp<DimensionedField<TypeR, GeoMesh, Field>> New
(
    const tmp<DimensionedField<TypeR, GeoMesh, Field>>& tdf1,
    const word& name,
    const dimensionSet& dimensions
)
{
    DimensionedField<TypeR, GeoMesh, Field>& df1 =
        const_cast<DimensionedField<TypeR, GeoMesh, Field>&>(tdf1());

    if (tdf1.isTmp())
    {
        df1.rename(name);
        df1.dimensions().reset(dimensions);
        return tdf1;
    }
    else
    {
        return DimensionedField<TypeR, GeoMesh, Field>::New
        (
            name,
            df1.mesh(),
            dimensions
        );
    }
}


template
<
    class TypeR,
    class Type1,
    class GeoMesh,
    template<class> class PrimitiveField
>
class reuseTmpDimensionedField
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh, Field>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh, PrimitiveField>>& tdf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const DimensionedField<Type1, GeoMesh, PrimitiveField>& df1 = tdf1();

        return DimensionedField<TypeR, GeoMesh, Field>::New
        (
            name,
            df1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class GeoMesh>
class reuseTmpDimensionedField<TypeR, TypeR, GeoMesh, Field>
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh, Field>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh, Field>>& tdf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        DimensionedField<TypeR, GeoMesh, Field>& df1 =
            const_cast<DimensionedField<TypeR, GeoMesh, Field>&>(tdf1());

        if (tdf1.isTmp())
        {
            df1.rename(name);
            df1.dimensions().reset(dimensions);
            return tdf1;
        }
        else
        {
            return DimensionedField<TypeR, GeoMesh, Field>::New
            (
                name,
                df1.mesh(),
                dimensions
            );
        }
    }
};


template
<
    class TypeR,
    class Type1,
    class Type2,
    class GeoMesh,
    template<class> class PrimitiveField1,
    template<class> class PrimitiveField2
>
class reuseTmpTmpDimensionedField
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh, Field>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh, PrimitiveField1>>& tdf1,
        const tmp<DimensionedField<Type2, GeoMesh, PrimitiveField2>>& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const DimensionedField<Type1, GeoMesh, PrimitiveField1>& df1 = tdf1();

        return DimensionedField<TypeR, GeoMesh, Field>::New
        (
            name,
            df1.mesh(),
            dimensions
        );
    }
};


template
<
    class TypeR,
    class Type1,
    class GeoMesh,
    template<class> class PrimitiveField1
>
class reuseTmpTmpDimensionedField
<TypeR, Type1, TypeR, GeoMesh, PrimitiveField1, Field>
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh, Field>> New
    (
        const tmp<DimensionedField<Type1, GeoMesh, PrimitiveField1>>& tdf1,
        const tmp<DimensionedField<TypeR, GeoMesh, Field>>& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        DimensionedField<TypeR, GeoMesh, Field>& df2 =
            const_cast<DimensionedField<TypeR, GeoMesh, Field>&>(tdf2());

        if (tdf2.isTmp())
        {
            df2.rename(name);
            df2.dimensions().reset(dimensions);
            return tdf2;
        }
        else
        {
            return DimensionedField<TypeR, GeoMesh, Field>::New
            (
                name,
                df2.mesh(),
                dimensions
            );
        }
    }
};


template
<
    class TypeR,
    class Type2,
    class GeoMesh,
    template<class> class PrimitiveField2
>
class reuseTmpTmpDimensionedField
<TypeR, TypeR, Type2, GeoMesh, Field, PrimitiveField2>
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh, Field>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh, Field>>& tdf1,
        const tmp<DimensionedField<Type2, GeoMesh, PrimitiveField2>>& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        DimensionedField<TypeR, GeoMesh, Field>& df1 =
            const_cast<DimensionedField<TypeR, GeoMesh, Field>&>(tdf1());

        if (tdf1.isTmp())
        {
            df1.rename(name);
            df1.dimensions().reset(dimensions);
            return tdf1;
        }
        else
        {
            return DimensionedField<TypeR, GeoMesh, Field>::New
            (
                name,
                df1.mesh(),
                dimensions
            );
        }
    }
};


template<class TypeR, class GeoMesh>
class reuseTmpTmpDimensionedField<TypeR, TypeR, TypeR, GeoMesh, Field, Field>
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh, Field>> New
    (
        const tmp<DimensionedField<TypeR, GeoMesh, Field>>& tdf1,
        const tmp<DimensionedField<TypeR, GeoMesh, Field>>& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        DimensionedField<TypeR, GeoMesh, Field>& df1 =
            const_cast<DimensionedField<TypeR, GeoMesh, Field>&>(tdf1());
        DimensionedField<TypeR, GeoMesh, Field>& df2 =
            const_cast<DimensionedField<TypeR, GeoMesh, Field>&>(tdf2());

        if (tdf1.isTmp())
        {
            df1.rename(name);
            df1.dimensions().reset(dimensions);
            return tdf1;
        }
        else if (tdf2.isTmp())
        {
            df2.rename(name);
            df2.dimensions().reset(dimensions);
            return tdf2;
        }
        else
        {
            return DimensionedField<TypeR, GeoMesh, Field>::New
            (
                name,
                df1.mesh(),
                dimensions
            );
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
