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

#ifndef DimensionedFieldReuseFunctions_H
#define DimensionedFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class TypeR, class Type1, class GeoMesh>
class reuseTmpDimensionedField
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh> > New
    (
        const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const DimensionedField<Type1, GeoMesh>& df1 = tdf1();

        return tmp<DimensionedField<TypeR, GeoMesh> >
        (
            new DimensionedField<TypeR, GeoMesh>
            (
                IOobject
                (
                    name,
                    df1.instance(),
                    df1.db()
                ),
                df1.mesh(),
                dimensions
            )
        );
    }

    static void clear(const tmp<DimensionedField<Type1, GeoMesh> >& tdf1)
    {
        tdf1.clear();
    }
};


template<class TypeR, class GeoMesh>
class reuseTmpDimensionedField<TypeR, TypeR, GeoMesh>
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh> > New
    (
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf1,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        DimensionedField<TypeR, GeoMesh>& df1 =
            const_cast<DimensionedField<TypeR, GeoMesh>& >(tdf1());

        if (tdf1.isTmp())
        {
            df1.rename(name);
            df1.dimensions().reset(dimensions);
            return tdf1;
        }
        else
        {
            return tmp<DimensionedField<TypeR, GeoMesh> >
            (
                new DimensionedField<TypeR, GeoMesh>
                (
                    IOobject
                    (
                        name,
                        df1.instance(),
                        df1.db()
                    ),
                    df1.mesh(),
                    dimensions
                )
            );
        }
    }

    static void clear(const tmp<DimensionedField<TypeR, GeoMesh> >& tdf1)
    {
        if (tdf1.isTmp())
        {
            tdf1.ptr();
        }
    }
};


template<class TypeR, class Type1, class Type12, class Type2, class GeoMesh>
class reuseTmpTmpDimensionedField
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh> > New
    (
        const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,
        const tmp<DimensionedField<Type2, GeoMesh> >& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const DimensionedField<Type1, GeoMesh>& df1 = tdf1();

        return tmp<DimensionedField<TypeR, GeoMesh> >
        (
            new DimensionedField<TypeR, GeoMesh>
            (
                IOobject
                (
                    name,
                    df1.instance(),
                    df1.db()
                ),
                df1.mesh(),
                dimensions
            )
        );
    }

    static void clear
    (
        const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,
        const tmp<DimensionedField<Type2, GeoMesh> >& tdf2
    )
    {
        tdf1.clear();
        tdf2.clear();
    }
};


template<class TypeR, class Type1, class Type12, class GeoMesh>
class reuseTmpTmpDimensionedField<TypeR, Type1, Type12, TypeR, GeoMesh>
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh> > New
    (
        const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        const DimensionedField<Type1, GeoMesh>& df1 = tdf1();
        DimensionedField<TypeR, GeoMesh>& df2 =
            const_cast<DimensionedField<TypeR, GeoMesh>& >(tdf2());

        if (tdf2.isTmp())
        {
            df2.rename(name);
            df2.dimensions().reset(dimensions);
            return tdf2;
        }
        else
        {
            return tmp<DimensionedField<TypeR, GeoMesh> >
            (
                new DimensionedField<TypeR, GeoMesh>
                (
                    IOobject
                    (
                        name,
                        df1.instance(),
                        df1.db()
                    ),
                    df1.mesh(),
                    dimensions
                )
            );
        }
    }

    static void clear
    (
        const tmp<DimensionedField<Type1, GeoMesh> >& tdf1,
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf2
    )
    {
        tdf1.clear();
        if (tdf2.isTmp())
        {
            tdf2.ptr();
        }
    }
};


template<class TypeR, class Type2, class GeoMesh>
class reuseTmpTmpDimensionedField<TypeR, TypeR, TypeR, Type2, GeoMesh>
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh> > New
    (
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf1,
        const tmp<DimensionedField<Type2, GeoMesh> >& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        DimensionedField<TypeR, GeoMesh>& df1 =
            const_cast<DimensionedField<TypeR, GeoMesh>& >(tdf1());

        if (tdf1.isTmp())
        {
            df1.rename(name);
            df1.dimensions().reset(dimensions);
            return tdf1;
        }
        else
        {
            return tmp<DimensionedField<TypeR, GeoMesh> >
            (
                new DimensionedField<TypeR, GeoMesh>
                (
                    IOobject
                    (
                        name,
                        df1.instance(),
                        df1.db()
                    ),
                    df1.mesh(),
                    dimensions
                )
            );
        }
    }

    static void clear
    (
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf1,
        const tmp<DimensionedField<Type2, GeoMesh> >& tdf2
    )
    {
        if (tdf1.isTmp())
        {
            tdf1.ptr();
        }
        tdf2.clear();
    }
};


template<class TypeR, class GeoMesh>
class reuseTmpTmpDimensionedField<TypeR, TypeR, TypeR, TypeR, GeoMesh>
{
public:

    static tmp<DimensionedField<TypeR, GeoMesh> > New
    (
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf1,
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf2,
        const word& name,
        const dimensionSet& dimensions
    )
    {
        DimensionedField<TypeR, GeoMesh>& df1 =
            const_cast<DimensionedField<TypeR, GeoMesh>& >(tdf1());
        DimensionedField<TypeR, GeoMesh>& df2 =
            const_cast<DimensionedField<TypeR, GeoMesh>& >(tdf2());

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
            return tmp<DimensionedField<TypeR, GeoMesh> >
            (
                new DimensionedField<TypeR, GeoMesh>
                (
                    IOobject
                    (
                        name,
                        df1.instance(),
                        df1.db()
                    ),
                    df1.mesh(),
                    dimensions
                )
            );
        }
    }

    static void clear
    (
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf1,
        const tmp<DimensionedField<TypeR, GeoMesh> >& tdf2
    )
    {
        if (tdf1.isTmp())
        {
            tdf1.ptr();
            tdf2.clear();
        }
        else if (tdf2.isTmp())
        {
            tdf1.clear();
            tdf2.ptr();
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
