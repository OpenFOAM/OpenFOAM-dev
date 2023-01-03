/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#ifndef GeometricFieldReuseFunctions_H
#define GeometricFieldReuseFunctions_H

#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh, template<class> class PrimitiveField>
bool reusable(const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>& tgf)
{
    typedef GeometricField<Type, GeoMesh, PrimitiveField> geoField;

    if (tgf.isTmp())
    {
        if (geoField::debug)
        {
            const geoField& gf = tgf();
            const typename geoField::Boundary& gbf = gf.boundaryField();

            forAll(gbf, patchi)
            {
                if
                (
                    !polyPatch::constraintType(gbf[patchi].patch().type())
                 && !isA<typename geoField::Patch::Calculated>(gbf[patchi])
                )
                {
                    WarningInFunction
                        << "Attempt to reuse temporary with non-reusable BC "
                        << gbf[patchi].type() << endl;

                    return false;
                }
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class TypeR, class GeoMesh>
tmp<GeometricField<TypeR, GeoMesh, Field>> New
(
    const tmp<GeometricField<TypeR, GeoMesh, Field>>& tgf1,
    const word& name,
    const dimensionSet& dimensions,
    const bool initRet = false
)
{
    GeometricField<TypeR, GeoMesh, Field>& gf1 =
        const_cast<GeometricField<TypeR, GeoMesh, Field>&>(tgf1());

    if (reusable(tgf1))
    {
        gf1.rename(name);
        gf1.dimensions().reset(dimensions);
        return tgf1;
    }
    else
    {
        tmp<GeometricField<TypeR, GeoMesh, Field>> rtgf
        (
            GeometricField<TypeR, GeoMesh, Field>::New
            (
                name,
                gf1.mesh(),
                dimensions
            )
        );

        if (initRet)
        {
            rtgf.ref() == tgf1();
        }

        return rtgf;
    }
}


template
<
    class TypeR,
    class Type1,
    class GeoMesh,
    template<class> class PrimitiveField
>
class reuseTmpGeometricField
{
public:

    static tmp<GeometricField<TypeR, GeoMesh, Field>> New
    (
        const tmp<GeometricField<Type1, GeoMesh, PrimitiveField>>& tgf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const GeometricField<Type1, GeoMesh, PrimitiveField>& gf1 = tgf1();

        return GeometricField<TypeR, GeoMesh, Field>::New
        (
            name,
            gf1.mesh(),
            dimensions
        );
    }
};


template<class TypeR, class GeoMesh>
class reuseTmpGeometricField<TypeR, TypeR, GeoMesh, Field>
{
public:

    static tmp<GeometricField<TypeR, GeoMesh, Field>> New
    (
        const tmp<GeometricField<TypeR, GeoMesh, Field>>& tgf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        GeometricField<TypeR, GeoMesh, Field>& gf1 =
            const_cast<GeometricField<TypeR, GeoMesh, Field>&>(tgf1());

        if (reusable(tgf1))
        {
            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            return tgf1;
        }
        else
        {
            return GeometricField<TypeR, GeoMesh, Field>::New
            (
                name,
                gf1.mesh(),
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
class reuseTmpTmpGeometricField
{
public:

    static tmp<GeometricField<TypeR, GeoMesh, Field>> New
    (
        const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tgf1,
        const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const GeometricField<Type1, GeoMesh, PrimitiveField1>& gf1 = tgf1();

        return GeometricField<TypeR, GeoMesh, PrimitiveField1>::New
        (
            name,
            gf1.mesh(),
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
class reuseTmpTmpGeometricField
<TypeR, Type1, TypeR, GeoMesh, PrimitiveField1, Field>
{
public:

    static tmp<GeometricField<TypeR, GeoMesh, Field>> New
    (
        const tmp<GeometricField<Type1, GeoMesh, PrimitiveField1>>& tgf1,
        const tmp<GeometricField<TypeR, GeoMesh, Field>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const GeometricField<Type1, GeoMesh, PrimitiveField1>& gf1 = tgf1();

        GeometricField<TypeR, GeoMesh, Field>& gf2 =
            const_cast<GeometricField<TypeR, GeoMesh, Field>&>(tgf2());

        if (reusable(tgf2))
        {
            gf2.rename(name);
            gf2.dimensions().reset(dimensions);
            return tgf2;
        }
        else
        {
            return GeometricField<TypeR, GeoMesh, Field>::New
            (
                name,
                gf1.mesh(),
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
class reuseTmpTmpGeometricField
<TypeR, TypeR, Type2, GeoMesh, Field, PrimitiveField2>
{
public:

    static tmp<GeometricField<TypeR, GeoMesh, Field>> New
    (
        const tmp<GeometricField<TypeR, GeoMesh, Field>>& tgf1,
        const tmp<GeometricField<Type2, GeoMesh, PrimitiveField2>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        GeometricField<TypeR, GeoMesh, Field>& gf1 =
            const_cast<GeometricField<TypeR, GeoMesh, Field>&>(tgf1());

        if (reusable(tgf1))
        {
            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            return tgf1;
        }
        else
        {
            return GeometricField<TypeR, GeoMesh, Field>::New
            (
                name,
                gf1.mesh(),
                dimensions
            );
        }
    }
};


template<class TypeR, class GeoMesh>
class reuseTmpTmpGeometricField<TypeR, TypeR, TypeR, GeoMesh, Field, Field>
{
public:

    static tmp<GeometricField<TypeR, GeoMesh, Field>> New
    (
        const tmp<GeometricField<TypeR, GeoMesh, Field>>& tgf1,
        const tmp<GeometricField<TypeR, GeoMesh, Field>>& tgf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        GeometricField<TypeR, GeoMesh, Field>& gf1 =
            const_cast<GeometricField<TypeR, GeoMesh, Field>&>(tgf1());
        GeometricField<TypeR, GeoMesh, Field>& gf2 =
            const_cast<GeometricField<TypeR, GeoMesh, Field>&>(tgf2());

        if (reusable(tgf1))
        {
            gf1.rename(name);
            gf1.dimensions().reset(dimensions);
            return tgf1;
        }
        else if (reusable(tgf2))
        {
            gf2.rename(name);
            gf2.dimensions().reset(dimensions);
            return tgf2;
        }
        else
        {
            return GeometricField<TypeR, GeoMesh, Field>::New
            (
                name,
                gf1.mesh(),
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
