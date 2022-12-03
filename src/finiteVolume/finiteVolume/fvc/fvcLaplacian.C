/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
tmp<VolField<Type>>
laplacian
(
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, scalar>::New
    (
        vf.mesh(),
        vf.mesh().schemes().laplacian(name)
    ).ref().fvcLaplacian(vf);
}


template<class Type>
tmp<VolField<Type>>
laplacian
(
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<VolField<Type>>
laplacian
(
    const VolField<Type>& vf
)
{
    return fvc::laplacian(vf, "laplacian(" + vf.name() + ')');
}


template<class Type>
tmp<VolField<Type>>
laplacian
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const VolField<Type>& vf,
    const word& name
)
{
    SurfaceField<GType> Gamma
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
tmp<VolField<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const VolField<Type>& vf
)
{
    SurfaceField<GType> Gamma
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
tmp<VolField<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().schemes().laplacian(name)
    ).ref().fvcLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf
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
tmp<VolField<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const VolField<Type>& vf
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
tmp<VolField<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const tmp<VolField<Type>>& tvf
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
tmp<VolField<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const tmp<VolField<Type>>& tvf
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
tmp<VolField<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::laplacianScheme<Type, GType>::New
    (
        vf.mesh(),
        vf.mesh().schemes().laplacian(name)
    ).ref().fvcLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>> laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const VolField<Type>& vf
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
tmp<VolField<Type>>
laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const VolField<Type>& vf
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(tgamma(), vf)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Laplacian
    (
        fvc::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolField<Type>> laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Laplacian
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
