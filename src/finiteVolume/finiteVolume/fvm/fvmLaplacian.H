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
    Calculate the matrix for the laplacian of the field.

SourceFiles
    fvmLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef fvmLaplacian_H
#define fvmLaplacian_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "zero.H"
#include "one.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvMatrix;

/*---------------------------------------------------------------------------*\
                      Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{
    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const VolField<Type>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const zero&,
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const zero&,
        const VolField<Type>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const one&,
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplacian
    (
        const one&,
        const VolField<Type>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&
    );


    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<fvMatrix<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> laplacianCorrection
    (
        const VolField<scalar>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplacianCorrection
    (
        const tmp<VolField<scalar>>&,
        const VolField<Type>&
    );


    template<class Type>
    tmp<fvMatrix<Type>> laplacianCorrection
    (
        const SurfaceField<scalar>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<fvMatrix<Type>> laplacianCorrection
    (
        const tmp<SurfaceField<scalar>>&,
        const VolField<Type>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvmLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
