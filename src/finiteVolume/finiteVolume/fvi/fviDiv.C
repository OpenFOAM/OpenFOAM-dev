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

#include "fviDiv.H"
#include "fvMesh.H"
#include "fviSurfaceIntegrate.H"
#include "divScheme.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvi
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<VolInternalField<Type>>
div
(
    const SurfaceField<Type>& ssf
)
{
    return VolInternalField<Type>::New
    (
        "div("+ssf.name()+')',
        fvi::surfaceIntegrate(ssf)
    );
}


template<class Type>
tmp<VolInternalField<Type>>
div
(
    const tmp<SurfaceField<Type>>& tssf
)
{
    tmp<VolInternalField<Type>> Div(fvi::div(tssf()));
    tssf.clear();
    return Div;
}


template<class Type>
tmp<VolInternalField<typename innerProduct<vector, Type>::type>>
div
(
    const VolField<Type>& vf,
    const word& name
)
{
    return fv::divScheme<Type>::New
    (
        vf.mesh(), vf.mesh().schemes().div(name)
    ).ref().fviDiv(vf);
}


template<class Type>
tmp<VolInternalField<typename innerProduct<vector, Type>::type>>
div
(
    const tmp<VolField<Type>>& tvvf,
    const word& name
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<VolInternalField<DivType>> Div
    (
        fvi::div(tvvf(), name)
    );
    tvvf.clear();
    return Div;
}

template<class Type>
tmp<VolInternalField<typename innerProduct<vector, Type>::type>>
div
(
    const VolField<Type>& vf
)
{
    return fvi::div(vf, "div("+vf.name()+')');
}


template<class Type>
tmp<VolInternalField<typename innerProduct<vector, Type>::type>>
div
(
    const tmp<VolField<Type>>& tvvf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<VolInternalField<DivType>> Div(fvi::div(tvvf()));
    tvvf.clear();
    return Div;
}


template<class Type>
tmp<VolInternalField<Type>>
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
    ).ref().fviDiv(flux, vf);
}


template<class Type>
tmp<VolInternalField<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Div
    (
        fvi::div(tflux(), vf, name)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<VolInternalField<Type>>
div
(
    const surfaceScalarField& flux,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Div
    (
        fvi::div(flux, tvf(), name)
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolInternalField<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const tmp<VolField<Type>>& tvf,
    const word& name
)
{
    tmp<VolInternalField<Type>> Div
    (
        fvi::div(tflux(), tvf(), name)
    );
    tflux.clear();
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolInternalField<Type>>
div
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf
)
{
    return fvi::div
    (
        flux, vf, "div("+flux.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<VolInternalField<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf
)
{
    tmp<VolInternalField<Type>> Div
    (
        fvi::div(tflux(), vf)
    );
    tflux.clear();
    return Div;
}


template<class Type>
tmp<VolInternalField<Type>>
div
(
    const surfaceScalarField& flux,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolInternalField<Type>> Div
    (
        fvi::div(flux, tvf())
    );
    tvf.clear();
    return Div;
}


template<class Type>
tmp<VolInternalField<Type>>
div
(
    const tmp<surfaceScalarField>& tflux,
    const tmp<VolField<Type>>& tvf
)
{
    tmp<VolInternalField<Type>> Div
    (
        fvi::div(tflux(), tvf())
    );
    tflux.clear();
    tvf.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvi

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
