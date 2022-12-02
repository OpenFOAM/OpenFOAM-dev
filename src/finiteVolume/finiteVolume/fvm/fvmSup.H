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

InNamespace
    Foam::fvm

Description
    Calculate the matrix for implicit and explicit sources.

SourceFiles
    fvmSup.C

\*---------------------------------------------------------------------------*/

#ifndef fvmSup_H
#define fvmSup_H

#include "volFieldsFwd.H"
#include "zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvMatrix;

/*---------------------------------------------------------------------------*\
                      Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{
   // Explicit source

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const DimensionedField<Type, volMesh>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const tmp<DimensionedField<Type, volMesh>>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const tmp<VolField<Type>>&,
            const VolField<Type>&
        );

        template<class Type>
        zeroField Su
        (
            const zero&,
            const VolField<Type>&
        );


    // Implicit source

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const volScalarField::Internal&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const tmp<volScalarField::Internal>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const tmp<volScalarField>&,
            const VolField<Type>&
        );


        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const dimensionedScalar&,
            const VolField<Type>&
        );


        template<class Type>
        zeroField Sp
        (
            const zero&,
            const VolField<Type>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const volScalarField::Internal&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const tmp<volScalarField::Internal>&,
            const VolField<Type>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const tmp<volScalarField>&,
            const VolField<Type>&
        );

        template<class Type>
        zeroField SuSp
        (
            const zero&,
            const VolField<Type>&
        );


    // Explicit and implicit sources provided as a pair

        template<class Type>
        tmp<fvMatrix<Type>> S
        (
            const Pair<tmp<volScalarField::Internal>>&,
            const VolField<Type>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvmSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
