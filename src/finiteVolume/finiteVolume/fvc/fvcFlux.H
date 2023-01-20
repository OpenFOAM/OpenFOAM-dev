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
    Calculate the face-flux of the given field.

SourceFiles
    fvcFluxTemplates.C

\*---------------------------------------------------------------------------*/


#ifndef fvcFlux_H
#define fvcFlux_H

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
    //- Return the face-flux field obtained from the given volVectorField
    template<class Type>
    tmp<SurfaceField<typename innerProduct<vector, Type>::type>> flux
    (
        const VolField<Type>& vf
    );

    //- Return the face-flux field obtained from the given tmp volVectorField
    template<class Type>
    tmp<SurfaceField<typename innerProduct<vector, Type>::type>>  flux
    (
        const tmp<VolField<Type>>& tvf
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const surfaceScalarField&,
        const tmp<VolField<Type>>&
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const surfaceScalarField&,
        const VolField<Type>&,
        Istream& schemeData
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const surfaceScalarField&,
        const VolField<Type>&,
        const word& name
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&,
        const word& name
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const surfaceScalarField&,
        const tmp<VolField<Type>>&,
        const word& name
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const tmp<surfaceScalarField>&,
        const tmp<VolField<Type>>&,
        const word& name
    );


    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const surfaceScalarField&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const tmp<surfaceScalarField>&,
        const VolField<Type>&
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const surfaceScalarField&,
        const tmp<VolField<Type>>&
    );

    template<class Type>
    tmp<SurfaceField<Type>> flux
    (
        const tmp<surfaceScalarField>&,
        const tmp<VolField<Type>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcFluxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
