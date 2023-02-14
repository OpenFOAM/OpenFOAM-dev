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

#include "MRFZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::MRFZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = this->Omega();

    const vectorField& Cfi = Cf;
    const vectorField& Sfi = Sf;
    scalarField& phii = phi.primitiveFieldRef();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];
        phii[facei] -= rho[facei]*(Omega ^ (Cfi[facei] - origin_)) & Sfi[facei];
    }

    makeRelativeRhoFlux(rho.boundaryField(), phi.boundaryFieldRef());
}


template<class RhoFieldType>
void Foam::MRFZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    FieldField<fvsPatchField, scalar>& phiBf
) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = this->Omega();

    forAll(patchFaces_, patchi)
    {
        if (!phiBf[patchi].fixesValue())
        {
            forAll(patchFaces_[patchi], i)
            {
                const label patchFacei = patchFaces_[patchi][i];

                phiBf[patchi][patchFacei] -=
                    rho[patchi][patchFacei]
                   *(Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin_))
                  & Sf.boundaryField()[patchi][patchFacei];
            }
        }
    }
}


template<class RhoFieldType>
void Foam::MRFZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    Field<scalar>& phi,
    const label patchi
) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = this->Omega();

    forAll(patchFaces_[patchi], i)
    {
        const label patchFacei = patchFaces_[patchi][i];

        phi[patchFacei] -=
            rho[patchFacei]
          * (Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin_))
          & Sf.boundaryField()[patchi][patchFacei];
    }
}


template<class RhoFieldType>
void Foam::MRFZone::makeAbsoluteRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = this->Omega();

    const vectorField& Cfi = Cf;
    const vectorField& Sfi = Sf;
    scalarField& phii = phi.primitiveFieldRef();

    // Internal faces
    forAll(internalFaces_, i)
    {
        const label facei = internalFaces_[i];
        phii[facei] += rho[facei]*(Omega ^ (Cfi[facei] - origin_)) & Sfi[facei];
    }

    surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

    forAll(patchFaces_, patchi)
    {
        if (!phiBf[patchi].fixesValue())
        {
            forAll(patchFaces_[patchi], i)
            {
                const label patchFacei = patchFaces_[patchi][i];

                phiBf[patchi][patchFacei] +=
                    rho.boundaryField()[patchi][patchFacei]
                   *(Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin_))
                  & Sf.boundaryField()[patchi][patchFacei];
            }
        }
    }
}


template<class Type>
void Foam::MRFZone::zero
(
    SurfaceField<Type>& phi
) const
{
    Field<Type>& phii = phi.primitiveFieldRef();

    forAll(internalFaces_, i)
    {
        phii[internalFaces_[i]] = Zero;
    }

    typename SurfaceField<Type>::Boundary& phibf =
        phi.boundaryFieldRef();

    forAll(patchFaces_, patchi)
    {
        forAll(patchFaces_[patchi], i)
        {
            phibf[patchi][patchFaces_[patchi][i]] = Zero;
        }
    }
}


// ************************************************************************* //
