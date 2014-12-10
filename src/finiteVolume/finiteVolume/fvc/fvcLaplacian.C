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

#include "fvcLaplacian.H"
#include "fvMesh.H"
#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, scalar>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().fvcLaplacian(vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::laplacian(vf, "laplacian(" + vf.name() + ')');
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    GeometricField<GType, fvsPatchField, surfaceMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fvc::laplacian(Gamma, vf, name);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const dimensioned<GType>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    GeometricField<GType, fvsPatchField, surfaceMesh> Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.instance(),
            vf.mesh(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return fvc::laplacian(Gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().fvcLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::laplacian
    (
        tgamma,
        vf,
        "laplacian(" + tgamma().name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvPatchField, volMesh>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    return fvc::laplacian
    (
        gamma,
        tvf,
        "laplacian(" + gamma.name() + ',' + tvf().name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvPatchField, volMesh> >& tgamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    return fvc::laplacian
    (
        tgamma,
        tvf,
        "laplacian(" + tgamma().name() + ',' + tvf().name() + ')'
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().laplacianScheme(name)
    )().fvcLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> > laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), vf)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> >
laplacian
(
    const GeometricField<GType, fvsPatchField, surfaceMesh>& gamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<GeometricField<Type, fvPatchField, volMesh> > laplacian
(
    const tmp<GeometricField<GType, fvsPatchField, surfaceMesh> >& tgamma,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > Laplacian
    (
        fvc::laplacian(tgamma(), tvf())
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
