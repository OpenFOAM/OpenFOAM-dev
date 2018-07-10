/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>>
Foam::fvc::scheme
(
    const surfaceScalarField& faceFlux,
    Istream& streamData
)
{
    return surfaceInterpolationScheme<Type>::New
    (
        faceFlux.mesh(),
        faceFlux,
        streamData
    );
}


template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>> Foam::fvc::scheme
(
    const surfaceScalarField& faceFlux,
    const word& name
)
{
    return surfaceInterpolationScheme<Type>::New
    (
        faceFlux.mesh(),
        faceFlux,
        faceFlux.mesh().interpolationScheme(name)
    );
}


template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>> Foam::fvc::scheme
(
    const fvMesh& mesh,
    Istream& streamData
)
{
    return surfaceInterpolationScheme<Type>::New
    (
        mesh,
        streamData
    );
}


template<class Type>
Foam::tmp<Foam::surfaceInterpolationScheme<Type>> Foam::fvc::scheme
(
    const fvMesh& mesh,
    const word& name
)
{
    return surfaceInterpolationScheme<Type>::New
    (
        mesh,
        mesh.interpolationScheme(name)
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, fvPatchField, volMesh> "
            << vf.name() << endl;
    }

    return scheme<Type>(faceFlux, schemeData)().interpolate(vf);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const surfaceScalarField& faceFlux,
    const word& name
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, fvPatchField, volMesh> "
            << vf.name() << " using " << name << endl;
    }

    return scheme<Type>(faceFlux, name)().interpolate(vf);
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const surfaceScalarField& faceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        interpolate(tvf(), faceFlux, name);

    tvf.clear();

    return tsf;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tFaceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        interpolate(vf, tFaceFlux(), name);

    tFaceFlux.clear();

    return tsf;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const tmp<surfaceScalarField>& tFaceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        interpolate(tvf(), tFaceFlux(), name);

    tvf.clear();
    tFaceFlux.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    Istream& schemeData
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, fvPatchField, volMesh> "
            << vf.name() << endl;
    }

    return scheme<Type>(vf.mesh(), schemeData)().interpolate(vf);
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, fvPatchField, volMesh> "
            << vf.name() << " using " << name
            << endl;
    }

    return scheme<Type>(vf.mesh(), name)().interpolate(vf);
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        interpolate(tvf(), name);

    tvf.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, fvPatchField, volMesh> "
            << vf.name() << " using run-time selected scheme"
            << endl;
    }

    return interpolate(vf, "interpolate(" + vf.name() + ')');
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvc::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsf =
        interpolate(tvf());
    tvf.clear();
    return tsf;
}


template<class Type>
Foam::tmp<Foam::FieldField<Foam::fvsPatchField, Type>>
Foam::fvc::interpolate
(
    const FieldField<fvPatchField, Type>& fvpff
)
{
    FieldField<fvsPatchField, Type>* fvspffPtr
    (
        new FieldField<fvsPatchField, Type>(fvpff.size())
    );

    forAll(*fvspffPtr, patchi)
    {
        fvspffPtr->set
        (
            patchi,
            fvsPatchField<Type>::NewCalculatedType(fvpff[patchi].patch()).ptr()
        );
        (*fvspffPtr)[patchi] = fvpff[patchi];
    }

    return tmp<FieldField<fvsPatchField, Type>>(fvspffPtr);
}


template<class Type>
Foam::tmp<Foam::FieldField<Foam::fvsPatchField, Type>>
Foam::fvc::interpolate
(
    const tmp<FieldField<fvPatchField, Type>>& tfvpff
)
{
    tmp<FieldField<fvsPatchField, Type>> tfvspff = interpolate(tfvpff());
    tfvpff.clear();
    return tfvspff;
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::innerProduct<Foam::vector, Type>::type,
        Foam::fvsPatchField,
        Foam::surfaceMesh
    >
>
Foam::fvc::dotInterpolate
(
    const surfaceVectorField& Sf,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, fvPatchField, volMesh> "
            << vf.name() << " using run-time selected scheme"
            << endl;
    }

    return scheme<Type>
    (
        vf.mesh(),
        "dotInterpolate(" + Sf.name() + ',' + vf.name() + ')'
    )().dotInterpolate(Sf, vf);
}


template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::innerProduct<Foam::vector, Type>::type,
        Foam::fvsPatchField,
        Foam::surfaceMesh
    >
>
Foam::fvc::dotInterpolate
(
    const surfaceVectorField& Sf,
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    tmp
    <
        GeometricField
        <
            typename Foam::innerProduct<Foam::vector, Type>::type,
            fvsPatchField,
            surfaceMesh
        >
    > tsf = dotInterpolate(Sf, tvf());
    tvf.clear();
    return tsf;
}


// ************************************************************************* //
