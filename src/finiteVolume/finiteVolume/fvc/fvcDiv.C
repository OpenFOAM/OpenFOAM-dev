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

#include "fvcDiv.H"
#include "fvMesh.H"
#include "fvcSurfaceIntegrate.H"
#include "divScheme.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolField<Type>>
div
(
    const SurfaceField<Type>& ssf
)
{
    return VolField<Type>::New
    (
        "div("+ssf.name()+')',
        fvc::surfaceIntegrate(ssf)
    );
}


template<class Type>
tmp<VolField<Type>>
div
(
    const tmp<SurfaceField<Type>>& tssf
)
{
    tmp<VolField<Type>> Div(fvc::div(tssf()));
    tssf.clear();
    return Div;
}


template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>>
div
(
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::divScheme<Type>::New
    (
        vf.mesh(), vf.mesh().schemes().div(name)
    ).ref().fvcDiv(vf);
}


template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>>
div
(
    const tmp<VolField<Type>>& tvvf,
    const word& name
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<VolField<DivType>> Div
    (
        fvc::div(tvvf(), name)
    );
    tvvf.clear();
    return Div;
}

template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>>
div
(
    const VolField<Type>& vf
)
{
    return fvc::div(vf, "div("+vf.name()+')');
}


template<class Type>
tmp<VolField<typename innerProduct<vector, Type>::type>>
div
(
    const tmp<VolField<Type>>& tvvf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<VolField<DivType>> Div(fvc::div(tvvf()));
    tvvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>>
div
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().schemes().div(name)
    ).ref().fvcDiv(flux, vf);
}


template<class Type>
tmp<VolField<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<VolField<Type>> Div
    (
        fvc::div(tflux(), vf, name)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>>
div
(
    const surfaceScalarField& flux,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Div
    (
        fvc::div(flux, tvf(), name)
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolField<Type>> Div
    (
        fvc::div(tflux(), tvf(), name)
    );
    tflux.clear();
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>>
div
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf
)
{
    return fvc::div
    (
        flux, vf, "div("+flux.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<VolField<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf
)
{
    tmp<VolField<Type>> Div
    (
        fvc::div(tflux(), vf)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>>
div
(
    const surfaceScalarField& flux,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Div
    (
        fvc::div(flux, tvf())
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolField<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolField<Type>> Div
    (
        fvc::div(tflux(), tvf())
    );
    tflux.clear();
    tvf.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
