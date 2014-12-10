/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<limitedSurfaceInterpolationScheme<Type> >
limitedSurfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (surfaceInterpolation::debug)
    {
        Info<< "limitedSurfaceInterpolationScheme<Type>::"
               "New(const fvMesh&, Istream&)"
               " : constructing limitedSurfaceInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "limitedSurfaceInterpolationScheme<Type>::"
            "New(const fvMesh&, Istream&)",
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
        FatalIOErrorIn
        (
            "limitedSurfaceInterpolationScheme<Type>::"
            "New(const fvMesh&, Istream&)",
            schemeData
        )   << "Unknown discretisation scheme "
            << schemeName << nl << nl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, schemeData);
}


// Return weighting factors for scheme given by name in dictionary
template<class Type>
tmp<limitedSurfaceInterpolationScheme<Type> >
limitedSurfaceInterpolationScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (surfaceInterpolation::debug)
    {
        Info<< "limitedSurfaceInterpolationScheme<Type>::New"
               "(const fvMesh&, const surfaceScalarField&, Istream&) : "
               "constructing limitedSurfaceInterpolationScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "limitedSurfaceInterpolationScheme<Type>::New"
            "(const fvMesh&, const surfaceScalarField&, Istream&)",
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
        FatalIOErrorIn
        (
            "limitedSurfaceInterpolationScheme<Type>::New"
            "(const fvMesh&, const surfaceScalarField&, Istream&)",
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
limitedSurfaceInterpolationScheme<Type>::~limitedSurfaceInterpolationScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<surfaceScalarField> limitedSurfaceInterpolationScheme<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    const surfaceScalarField& CDweights,
    tmp<surfaceScalarField> tLimiter
) const
{
    // Note that here the weights field is initialised as the limiter
    // from which the weight is calculated using the limiter value
    surfaceScalarField& Weights = tLimiter();

    scalarField& pWeights = Weights.internalField();

    forAll(pWeights, face)
    {
        pWeights[face] =
            pWeights[face]*CDweights[face]
          + (1.0 - pWeights[face])*pos(faceFlux_[face]);
    }

    surfaceScalarField::GeometricBoundaryField& bWeights =
        Weights.boundaryField();

    forAll(bWeights, patchI)
    {
        scalarField& pWeights = bWeights[patchI];

        const scalarField& pCDweights = CDweights.boundaryField()[patchI];
        const scalarField& pFaceFlux = faceFlux_.boundaryField()[patchI];

        forAll(pWeights, face)
        {
            pWeights[face] =
                pWeights[face]*pCDweights[face]
              + (1.0 - pWeights[face])*pos(pFaceFlux[face]);
        }
    }

    return tLimiter;
}

template<class Type>
tmp<surfaceScalarField> limitedSurfaceInterpolationScheme<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
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
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
limitedSurfaceInterpolationScheme<Type>::flux
(
    const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
    return faceFlux_*this->interpolate(phi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
