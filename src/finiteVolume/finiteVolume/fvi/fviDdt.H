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

InNamespace
    Foam::fvi

Description
    Calculate the first temporal derivative returning a volInternalField

SourceFiles
    fviDdt.C

\*---------------------------------------------------------------------------*/


#ifndef fviDdt_H
#define fviDdt_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedTypes.H"
#include "one.H"
#include "geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvi functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvi
{
    template<class Type>
    tmp<VolInternalField<Type>> ddt(const dimensioned<Type>, const fvMesh&);

    template<class Type>
    tmp<VolInternalField<Type>> ddt(const VolInternalField<Type>&);

    template<class Type>
    tmp<VolInternalField<Type>> ddt(const VolField<Type>&);

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const dimensionedScalar&,
        const VolInternalField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const dimensionedScalar&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const VolInternalField<scalar>&,
        const VolInternalField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const volScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const one&,
        const VolInternalField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const one&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const VolInternalField<scalar>&,
        const VolInternalField<scalar>&,
        const VolInternalField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const volScalarField&,
        const volScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const one&,
        const one&,
        const VolInternalField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const one&,
        const one&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const one&,
        const VolInternalField<scalar>&,
        const VolInternalField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const one&,
        const volScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const VolInternalField<scalar>&,
        const one&,
        const VolInternalField<Type>&
    );

    template<class Type>
    tmp<VolInternalField<Type>> ddt
    (
        const volScalarField&,
        const one&,
        const VolField<Type>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
   #include "fviDdt.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
