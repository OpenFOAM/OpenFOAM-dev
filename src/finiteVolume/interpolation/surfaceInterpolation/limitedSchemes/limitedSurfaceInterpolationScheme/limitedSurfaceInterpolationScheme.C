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

#include "limitedSurfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvPatchField.H"

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::limitedSurfaceInterpolationScheme<Type>>
Foam::limitedSurfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Constructing limitedSurfaceInterpolationScheme<Type>" << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    typename MeshConstructorTable::iterator constructorIter =
        MeshConstructorTablePtr_->find(schemeName);

    if (constructorIter == MeshConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, schemeData);
}


template<class Type>
Foam::tmp<Foam::limitedSurfaceInterpolationScheme<Type>>
Foam::limitedSurfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (surfaceInterpolation::debug)
    {
        InfoInFunction
            << "Constructing limitedSurfaceInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    typename MeshFluxConstructorTable::iterator constructorIter =
        MeshFluxConstructorTablePtr_->find(schemeName);

    if (constructorIter == MeshFluxConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction
        (
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << "Valid schemes are :" << endl
            << MeshFluxConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::limitedSurfaceInterpolationScheme<Type>::
~limitedSurfaceInterpolationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::limitedSurfaceInterpolationScheme<Type>::weights
(
    const VolField<Type>& phi,
    const surfaceScalarField& CDweights,
    tmp<surfaceScalarField> tLimiter
) const
{
    // Note that here the weights field is initialised as the limiter
    // from which the weight is calculated using the limiter value
    surfaceScalarField& Weights = tLimiter.ref();

    scalarField& pWeights = Weights.primitiveFieldRef();

    forAll(pWeights, face)
    {
        pWeights[face] =
            pWeights[face]*CDweights[face]
          + (1.0 - pWeights[face])*pos0(faceFlux_[face]);
    }

    surfaceScalarField::Boundary& bWeights =
        Weights.boundaryFieldRef();

    forAll(bWeights, patchi)
    {
        scalarField& pWeights = bWeights[patchi];

        const scalarField& pCDweights = CDweights.boundaryField()[patchi];
        const scalarField& pFaceFlux = faceFlux_.boundaryField()[patchi];

        forAll(pWeights, face)
        {
            pWeights[face] =
                pWeights[face]*pCDweights[face]
              + (1.0 - pWeights[face])*pos0(pFaceFlux[face]);
        }
    }

    return tLimiter;
}

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::limitedSurfaceInterpolationScheme<Type>::weights
(
    const VolField<Type>& phi
) const
{
    return this->weights
    (
        phi,
        this->mesh().surfaceInterpolation::weights(),
        this->limiter(phi)
    );
}

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::limitedSurfaceInterpolationScheme<Type>::flux
(
    const VolField<Type>& phi
) const
{
    return faceFlux_*this->interpolate(phi);
}


// ************************************************************************* //
