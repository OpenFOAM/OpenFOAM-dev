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

InNamespace
    Foam::fvc

Description
    Calculate the first temporal derivative.

SourceFiles
    fvcDdt.C

\*---------------------------------------------------------------------------*/


#ifndef fvcDdt_H
#define fvcDdt_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedTypes.H"
#include "one.H"
#include "geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const dimensioned<Type>,
        const fvMesh&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const dimensionedScalar&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const volScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const one&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const volScalarField&,
        const volScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const one&,
        const one&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const one&,
        const volScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> ddt
    (
        const volScalarField&,
        const one&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<SurfaceField<Type>> ddt
    (
        const SurfaceField<Type>&
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const VolField<Type>& U,
        const SurfaceField<Type>& Uf
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi,
        const autoPtr<SurfaceField<Type>>& Uf
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<Type>& rhoUf
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi,
        const autoPtr<SurfaceField<Type>>& rhoUf
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<Type>& Uf
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi
    );

    template<class Type>
    tmp<SurfaceField<typename Foam::flux<Type>::type>> ddtCorr
    (
        const volScalarField& alpha,
        const volScalarField& rho,
        const VolField<Type>& U,
        const SurfaceField<typename Foam::flux<Type>::type>& phi,
        const autoPtr<SurfaceField<Type>>& Uf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcDdt.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
