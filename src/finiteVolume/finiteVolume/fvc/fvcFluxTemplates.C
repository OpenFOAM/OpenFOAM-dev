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

#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<SurfaceField<typename innerProduct<vector, Type>::type>> flux
(
    const VolField<Type>& vf
)
{
    return scheme<Type>
    (
        vf.mesh(),
        "flux(" + vf.name() + ')'
    )().dotInterpolate(vf.mesh().Sf(), vf);
}


template<class Type>
tmp<SurfaceField<typename innerProduct<vector, Type>::type>>  flux
(
    const tmp<VolField<Type>>& tvf
)
{
    tmp<SurfaceField<typename innerProduct<vector, Type>::type>> Flux
    (
        fvc::flux(tvf())
    );
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const VolField<Type>& vf,
    Istream& schemeData
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        phi,
        schemeData
    )().flux(phi, vf);
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const VolField<Type>& vf,
    const word& name
)
{
    return fvc::flux(phi, vf, vf.mesh().schemes().div(name));
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(tphi(), vf, name)
    );
    tphi.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(phi, tvf(), name)
    );
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(tphi(), tvf(), name)
    );
    tphi.clear();
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const VolField<Type>& vf
)
{
    return fvc::flux
    (
        phi, vf, "flux("+phi.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const VolField<Type>& vf
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(tphi(), vf)
    );
    tphi.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const surfaceScalarField& phi,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(phi, tvf())
    );
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<SurfaceField<Type>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<SurfaceField<Type>> Flux
    (
        fvc::flux(tphi(), tvf())
    );
    tphi.clear();
    tvf.clear();
    return Flux;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
