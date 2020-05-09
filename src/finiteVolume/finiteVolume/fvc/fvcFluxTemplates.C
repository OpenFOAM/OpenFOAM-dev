/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "fvMesh.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
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
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fvc::flux(phi, vf, vf.mesh().divScheme(name));
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> Flux
    (
        fvc::flux(tphi(), vf, name)
    );
    tphi.clear();
    return Flux;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const surfaceScalarField& phi,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> Flux
    (
        fvc::flux(phi, tvf(), name)
    );
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> Flux
    (
        fvc::flux(tphi(), tvf(), name)
    );
    tphi.clear();
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const surfaceScalarField& phi,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::flux
    (
        phi, vf, "flux("+phi.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> Flux
    (
        fvc::flux(tphi(), vf)
    );
    tphi.clear();
    return Flux;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const surfaceScalarField& phi,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> Flux
    (
        fvc::flux(phi, tvf())
    );
    tvf.clear();
    return Flux;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
flux
(
    const tmp<surfaceScalarField>& tphi,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> Flux
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
