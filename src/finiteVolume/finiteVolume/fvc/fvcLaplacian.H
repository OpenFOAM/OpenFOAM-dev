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
    Calculate the laplacian of the given field.

SourceFiles
    fvcLaplacian.C

\*---------------------------------------------------------------------------*/


#ifndef fvcLaplacian_H
#define fvcLaplacian_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvc functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvc
{
    template<class Type>
    tmp<VolField<Type>> laplacian
    (
        const VolField<Type>&,
        const word&
    );

    template<class Type>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type>
    tmp<VolField<Type>> laplacian
    (
        const VolField<Type>&
    );

    template<class Type>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<Type>>&
    );


    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const dimensioned<GType>&,
        const tmp<VolField<Type>>&
    );


    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const VolField<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const VolField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const VolField<GType>&,
        const tmp<VolField<Type>>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<VolField<GType>>&,
        const tmp<VolField<Type>>&
    );


    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const tmp<VolField<Type>>&,
        const word&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const VolField<Type>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const SurfaceField<GType>&,
        const tmp<VolField<Type>>&
    );

    template<class Type, class GType>
    tmp<VolField<Type>> laplacian
    (
        const tmp<SurfaceField<GType>>&,
        const tmp<VolField<Type>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
