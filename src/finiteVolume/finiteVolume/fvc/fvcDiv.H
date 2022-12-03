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
    Foam::fvc

Description
    Calculate the divergence of the given field.

SourceFiles
    fvcDiv.C

\*---------------------------------------------------------------------------*/


#ifndef fvcDiv_H
#define fvcDiv_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<VolField<Type>> div
    (
        const SurfaceField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> div
    (
        const tmp<SurfaceField<Type>>&
    );


    template<class Type>
    tmp<VolField<typename innerProduct<vector, Type>::type>> div
    (
        const VolField<Type>&,
        const word& name
    );

    template<class Type>
    tmp<VolField<typename innerProduct<vector, Type>::type>> div
    (
        const tmp<VolField<Type>>&,
        const word& name
    );


    template<class Type>
    tmp<VolField<typename innerProduct<vector, Type>::type>> div
    (
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<typename innerProduct<vector, Type>::type>> div
    (
        const tmp<VolField<Type>>&
    );


    template<class Type>
    tmp<VolField<Type>> div
    (
        const surfaceScalarField&,
        const VolField<Type>&,
        const word& name
    );

    template<class Type>
    tmp<VolField<Type>> div
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&,
        const word& name
    );

    template<class Type>
    tmp<VolField<Type>> div
    (
        const surfaceScalarField&,
        const tmp<VolField<Type>>&,
        const word& name
    );

    template<class Type>
    tmp<VolField<Type>> div
    (
        const tmp<surfaceScalarField>&,
        const tmp<VolField<Type>>&,
        const word& name
    );


    template<class Type>
    tmp<VolField<Type>> div
    (
        const surfaceScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> div
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> div
    (
        const surfaceScalarField&,
        const tmp<VolField<Type>>&
    );

    template<class Type>
    tmp<VolField<Type>> div
    (
        const tmp<surfaceScalarField>&,
        const tmp<VolField<Type>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcDiv.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
