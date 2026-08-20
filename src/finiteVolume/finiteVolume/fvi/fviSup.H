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
    Calculate the field for explicit evaluation of implicit and explicit
    sources returning a volInternalField

SourceFiles
    fviSup.C

\*---------------------------------------------------------------------------*/

#ifndef fviSup_H
#define fviSup_H

#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvi functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvi
{
   // Explicit source

        template<class Type>
        tmp<VolInternalField<Type>> Su
        (
            const VolInternalField<Type>&,
            const VolField<Type>&
        );

        template
        <
            class Type,
            class Expression,
            class = EnableIfExpressionable<Expression>
        >
        tmp<VolInternalField<Type>> Su
        (
            const Expression& e,
            const VolField<Type>&
        );

        template<class Type>
        tmp<VolInternalField<Type>> Su
        (
            const tmp<VolInternalField<Type>>&,
            const VolField<Type>&
        );


    // Implicit source

        template<class Type>
        tmp<VolInternalField<Type>> Sp
        (
            const VolInternalField<scalar>&,
            const VolField<Type>&
        );

        template
        <
            class Type,
            class Expression,
            class = EnableIfExpressionable<Expression>
        >
        tmp<VolInternalField<Type>> Sp
        (
            const Expression& e,
            const VolField<Type>&
        );

        template<class Type>
        tmp<VolInternalField<Type>> Sp
        (
            const tmp<VolInternalField<scalar>>&,
            const VolField<Type>&
        );


        template<class Type>
        tmp<VolInternalField<Type>> Sp
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        template<class Type>
        tmp<VolInternalField<Type>> SuSp
        (
            const VolInternalField<scalar>&,
            const VolField<Type>&
        );

        template
        <
            class Type,
            class Expression,
            class = EnableIfExpressionable<Expression>
        >
        tmp<VolInternalField<Type>> SuSp
        (
            const Expression& e,
            const VolField<Type>&
        );

        template<class Type>
        tmp<VolInternalField<Type>> SuSp
        (
            const tmp<VolInternalField<scalar>>&,
            const VolField<Type>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fviSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
