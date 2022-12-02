/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

#include "boundedConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<SurfaceField<Type>>
boundedConvectionScheme<Type>::interpolate
(
    const surfaceScalarField& phi,
    const VolField<Type>& vf
) const
{
    return scheme_().interpolate(phi, vf);
}


template<class Type>
tmp<SurfaceField<Type>>
boundedConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const VolField<Type>& vf
) const
{
    return scheme_().flux(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
boundedConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const VolField<Type>& vf
) const
{
    return
        scheme_().fvmDiv(faceFlux, vf)
      - fvm::Sp(fvc::surfaceIntegrate(faceFlux), vf);
}


template<class Type>
tmp<VolField<Type>>
boundedConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const VolField<Type>& vf
) const
{
    return
        scheme_().fvcDiv(faceFlux, vf)
      - fvc::surfaceIntegrate(faceFlux)*vf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
