/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "phaseStabilised.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::volScalarField& Foam::phaseStabilised<Type>::lookupAlpha
(
    const surfaceScalarField& faceFlux
) const
{
    return faceFlux.db().lookupObject<volScalarField>
    (
        IOobject::groupName("alpha", faceFlux.group())
    );
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::phaseStabilised<Type>::lambdaf() const
{
    return pos(upwind_.interpolate(alpha_) - 1e-3);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class Type>
Foam::phaseStabilised<Type>::phaseStabilised
(
    const fvMesh& mesh,
    Istream& is
)
:
    surfaceInterpolationScheme<Type>(mesh),
    faceFlux_
    (
        mesh.lookupObject<surfaceScalarField>
        (
            word(is)
        )
    ),
    tScheme_
    (
        surfaceInterpolationScheme<Type>::New(mesh, faceFlux_, is)
    ),
    upwind_(mesh, faceFlux_),
    alpha_(lookupAlpha(faceFlux_))
{}


template<class Type>
Foam::phaseStabilised<Type>::phaseStabilised
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& is
)
:
    surfaceInterpolationScheme<Type>(mesh),
    faceFlux_(faceFlux),
    tScheme_
    (
        surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
    ),
    upwind_(mesh, faceFlux_),
    alpha_(lookupAlpha(faceFlux_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::phaseStabilised<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const surfaceScalarField lambdaf(this->lambdaf());
    return lambdaf*tScheme_().weights(vf) + (1 - lambdaf)*upwind_.weights();
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::phaseStabilised<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (tScheme_().corrected())
    {
        return lambdaf()*tScheme_().correction(vf);
    }
    else
    {
        return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        (
            nullptr
        );
    }
}


// ************************************************************************* //
