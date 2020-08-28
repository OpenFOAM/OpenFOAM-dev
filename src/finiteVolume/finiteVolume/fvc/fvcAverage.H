/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Area-weighted average a surfaceField creating a volField

SourceFiles
    fvcAverage.C

\*---------------------------------------------------------------------------*/


#ifndef fvcAverage_H
#define fvcAverage_H

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
    //- Area-weighted average a surfaceField creating a volField
    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh>> average
    (
        const GeometricField<Type, fvsPatchField, surfaceMesh>&
    );

    //- Area-weighted average a surfaceField creating a volField
    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh>> average
    (
        const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&
    );


    //- Area-weighted average a surfaceField creating a volField
    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh>> average
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    //- Area-weighted average a surfaceField creating a volField
    template<class Type>
    tmp<GeometricField<Type, fvPatchField, volMesh>> average
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvcAverage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
