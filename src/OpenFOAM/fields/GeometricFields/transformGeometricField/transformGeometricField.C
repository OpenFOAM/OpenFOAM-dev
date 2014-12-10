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

Description
    Spatial transformation functions for FieldFields.

\*---------------------------------------------------------------------------*/

#include "transformGeometricField.H"
#include "transformField.H"
#include "transformFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
    GeometricField<Type, PatchField, GeoMesh>& rtf,
    const GeometricField<tensor, PatchField, GeoMesh>& trf,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    transform(rtf.internalField(), trf.internalField(), tf.internalField());
    transform(rtf.boundaryField(), trf.boundaryField(), tf.boundaryField());
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
    const GeometricField<tensor, PatchField, GeoMesh>& trf,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > tranf
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                "transform(" + trf.name() + ',' + tf.name() + ')',
                tf.instance(),
                tf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tf.mesh(),
            tf.dimensions()
        )
    );

    transform(tranf(), trf, tf);

    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
    const GeometricField<tensor, PatchField, GeoMesh>& trf,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& ttf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > tranf =
        transform(trf, ttf());
    ttf.clear();
    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& ttrf,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > tranf =
        transform(ttrf(), tf);
    ttrf.clear();
    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
    const tmp<GeometricField<tensor, PatchField, GeoMesh> >& ttrf,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& ttf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > tranf =
        transform(ttrf(), ttf());
    ttf.clear();
    ttrf.clear();
    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void transform
(
    GeometricField<Type, PatchField, GeoMesh>& rtf,
    const dimensionedTensor& t,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    transform(rtf.internalField(), t.value(), tf.internalField());
    transform(rtf.boundaryField(), t.value(), tf.boundaryField());
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
    const dimensionedTensor& t,
    const GeometricField<Type, PatchField, GeoMesh>& tf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > tranf
    (
        new GeometricField<vector, PatchField, GeoMesh>
        (
            IOobject
            (
                "transform(" + t.name() + ',' + tf.name() + ')',
                tf.instance(),
                tf.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            tf.mesh(),
            tf.dimensions()
        )
    );

    transform(tranf(), t, tf);

    return tranf;
}


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh> > transform
(
    const dimensionedTensor& t,
    const tmp<GeometricField<Type, PatchField, GeoMesh> >& ttf
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh> > tranf =
        transform(t, ttf());
    ttf.clear();
    return tranf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
