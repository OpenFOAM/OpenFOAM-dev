/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "fvmDiv.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "convectionScheme.H"
#include "fvmSup.H"
#include "fviDiv.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::div
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
    )().fvmDiv(flux, vf);
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf,
    const word& name
)
{
    tmp<fvMatrix<Type>> tdiv(fvm::div(tflux(), vf, name));
    tflux.clear();
    return tdiv;
}


template<class Type, class Expression, class>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::div
(
    const Expression& e,
    const VolField<Type>& vf,
    const word& name
)
{
    return div(eval(e), vf, name);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::div
(
    const surfaceScalarField& flux,
    const VolField<Type>& vf
)
{
    return fvm::div(flux, vf, "div("+flux.name()+','+vf.name()+')');
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::div
(
    const tmp<surfaceScalarField>& tflux,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tdiv(fvm::div(tflux(), vf));
    tflux.clear();
    return tdiv;
}


template<class Type, class Expression, class>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::div
(
    const Expression& e,
    const VolField<Type>& vf
)
{
    return div(eval(e), vf);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::divc
(
    const tmp<SurfaceField<Type>>& tflux,
    const VolField<Type>& vf
)
{
    tmp<fvMatrix<Type>> tdivc(fvm::Su(fvi::div(tflux()), vf));

    if (vf.mesh().schemes().fluxRequired(vf.name()))
    {
        tdivc.ref().faceFluxCorrectionPtr() = tflux.ptr();
    }
    else
    {
        tflux.clear();
    }

    return tdivc;
}


template<class Type, class Expression, class>
Foam::tmp<Foam::fvMatrix<Type>>
Foam::fvm::divc
(
    const Expression& e,
    const VolField<Type>& vf
)
{
    return divc(eval(e), vf);
}


// ************************************************************************* //
