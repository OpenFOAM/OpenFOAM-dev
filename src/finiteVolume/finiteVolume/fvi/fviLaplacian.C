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

#include "fviLaplacian.H"
#include "fvMesh.H"
#include "laplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvi
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolInternalField<Type>>
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
    ).ref().fviLaplacian(vf);
}


template<class Type>
tmp<VolInternalField<Type>>
laplacian
(
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type>
tmp<VolInternalField<Type>>
laplacian
(
    const VolField<Type>& vf
)
{
    return fvi::laplacian(vf, "laplacian(" + vf.name() + ')');
}


template<class Type>
tmp<VolInternalField<Type>>
laplacian
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<VolInternalField<Type>>
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

    return fvi::laplacian(Gamma, vf, name);
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
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

    return fvi::laplacian(Gamma, vf);
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const dimensioned<GType>& gamma,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<VolInternalField<Type>>
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
    ).ref().fviLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvi::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const VolField<Type>& vf
)
{
    return fvi::laplacian
    (
        tgamma,
        vf,
        "laplacian(" + tgamma().name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const VolField<GType>& gamma,
    const tmp<VolField<Type>>& tvf
)
{
    return fvi::laplacian
    (
        gamma,
        tvf,
        "laplacian(" + gamma.name() + ',' + tvf().name() + ')'
    );
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const tmp<VolField<GType>>& tgamma,
    const tmp<VolField<Type>>& tvf
)
{
    return fvi::laplacian
    (
        tgamma,
        tvf,
        "laplacian(" + tgamma().name() + ',' + tvf().name() + ')'
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class GType>
tmp<VolInternalField<Type>>
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
    ).ref().fviLaplacian(gamma, vf);
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(tgamma(), vf, name)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(gamma, tvf(), name)
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>> laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(tgamma(), tvf(), name)
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const VolField<Type>& vf
)
{
    return fvi::laplacian
    (
        gamma,
        vf,
        "laplacian(" + gamma.name() + ',' + vf.name() + ')'
    );
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const VolField<Type>& vf
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(tgamma(), vf)
    );
    tgamma.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>>
laplacian
(
    const SurfaceField<GType>& gamma,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(gamma, tvf())
    );
    tvf.clear();
    return Laplacian;
}


template<class Type, class GType>
tmp<VolInternalField<Type>> laplacian
(
    const tmp<SurfaceField<GType>>& tgamma,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolInternalField<Type>> Laplacian
    (
        fvi::laplacian(tgamma(), tvf())
    );
    tgamma.clear();
    tvf.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
