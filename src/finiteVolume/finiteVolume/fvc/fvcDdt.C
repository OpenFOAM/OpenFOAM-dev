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

#include "fvcDdt.H"
#include "fvMesh.H"
#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>>
ddt
(
    const dimensioned<Type> dt,
    const fvMesh& mesh
)
{
    return fv::ddtScheme<Type>::New
    (
        mesh,
        mesh.schemes().ddt("ddt(" + dt.name() + ')')
    ).ref().fvcDdt(dt);
}


template<class Type>
tmp<VolField<Type>>
ddt
(
    const VolField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + vf.name() + ')')
    ).ref().fvcDdt(vf);
}


template<class Type>
tmp<VolField<Type>>
ddt
(
    const dimensionedScalar& rho,
    const VolField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<VolField<Type>>
ddt
(
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().fvcDdt(rho, vf);
}


template<class Type>
tmp<VolField<Type>>
ddt
(
    const one&,
    const VolField<Type>& vf
)
{
    return ddt(vf);
}


template<class Type>
tmp<VolField<Type>>
ddt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return fv::ddtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemes().ddt
        (
            "ddt("
          + alpha.name() + ','
          + rho.name() + ','
          + vf.name() + ')'
        )
    ).ref().fvcDdt(alpha, rho, vf);
}


template<class Type>
tmp<VolField<Type>>
ddt
(
    const one&,
    const one&,
    const VolField<Type>& vf
)
{
    return ddt(vf);
}


template<class Type>
tmp<VolField<Type>>
ddt
(
    const one&,
    const volScalarField& rho,
    const VolField<Type>& vf
)
{
    return ddt(rho, vf);
}


template<class Type>
tmp<VolField<Type>>
ddt
(
    const volScalarField& alpha,
    const one&,
    const VolField<Type>& vf
)
{
    return ddt(alpha, vf);
}


template<class Type>
tmp<SurfaceField<Type>>
ddt
(
    const SurfaceField<Type>& sf
)
{
    return fv::ddtScheme<Type>::New
    (
        sf.mesh(),
        sf.mesh().schemes().ddt("ddt(" + sf.name() + ')')
    ).ref().fvcDdt(sf);
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf
)
{
    return fv::ddtScheme<Type>::New
    (
        U.mesh(),
        U.mesh().schemes().ddt("ddt(" + U.name() + ')')
    ).ref().fvcDdtUfCorr(U, Uf);
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const VolField<Type>& U,
    const SurfaceField<typename Foam::flux<Type>::type>& phi
)
{
    return fv::ddtScheme<Type>::New
    (
        U.mesh(),
        U.mesh().schemes().ddt("ddt(" + U.name() + ')')
    ).ref().fvcDdtPhiCorr(U, phi);
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const VolField<Type>& U,
    const SurfaceField<typename Foam::flux<Type>::type>& phi,
    const autoPtr<SurfaceField<Type>>& Uf
)
{
    if (Uf.valid())
    {
        return ddtCorr(U, Uf());
    }
    else
    {
        return ddtCorr(U, phi);
    }
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<Type>& rhoUf
)
{
    return fv::ddtScheme<Type>::New
    (
        U.mesh(),
        U.mesh().schemes().ddt
        (
            "ddt(" + rho.name() + U.name() + ')'
        )
    ).ref().fvcDdtUfCorr(rho, U, rhoUf);
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<typename Foam::flux<Type>::type>& phi
)
{
    return fv::ddtScheme<Type>::New
    (
        U.mesh(),
        U.mesh().schemes().ddt("ddt(" + rho.name() + ',' + U.name() + ')')
    ).ref().fvcDdtPhiCorr(rho, U, phi);
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<typename Foam::flux<Type>::type>& phi,
    const autoPtr<SurfaceField<Type>>& rhoUf
)
{
    if (rhoUf.valid())
    {
        return ddtCorr(rho, U, rhoUf());
    }
    else
    {
        return ddtCorr(rho, U, phi);
    }
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<Type>& Uf
)
{
    return fv::ddtScheme<Type>::New
    (
        U.mesh(),
        U.mesh().schemes().ddt
        (
            "ddt(" + alpha.name() + rho.name() + ',' + U.name() + ')'
        )
    ).ref().fvcDdtUfCorr(alpha, rho, U, Uf);
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<typename Foam::flux<Type>::type>& phi
)
{
    return fv::ddtScheme<Type>::New
    (
        U.mesh(),
        U.mesh().schemes().ddt
        (
            "ddt(" + alpha.name() + rho.name() + ',' + U.name() + ')'
        )
    ).ref().fvcDdtPhiCorr(alpha, rho, U, phi);
}


template<class Type>
tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const VolField<Type>& U,
    const SurfaceField<typename Foam::flux<Type>::type>& phi,
    const autoPtr<SurfaceField<Type>>& Uf
)
{
    if (Uf.valid())
    {
        return ddtCorr(alpha, rho, U, Uf());
    }
    else
    {
        return ddtCorr(alpha, rho, U, phi);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
