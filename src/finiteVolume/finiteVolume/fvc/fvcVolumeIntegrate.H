/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Volume integrate volField creating a volField.

    Volume integrate volField over the whole domain
    creating a dimensioned\<Type\>

SourceFiles
    fvcVolumeIntegrate.C

\*---------------------------------------------------------------------------*/


#ifndef fvcVolumeIntegrate_H
#define fvcVolumeIntegrate_H

#include "volFieldsFwd.H"
#include "primitiveFieldsFwd.H"
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
    tmp<Field<Type> > volumeIntegrate
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    tmp<Field<Type> > volumeIntegrate
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh> >&
    );


    template<class Type>
    tmp<Field<Type> > volumeIntegrate
    (
        const DimensionedField<Type, volMesh>&
    );

    template<class Type>
    tmp<Field<Type> > volumeIntegrate
    (
        const tmp<DimensionedField<Type, volMesh> >&
    );


    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const GeometricField<Type, fvPatchField, volMesh>&
    );

    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const tmp<GeometricField<Type, fvPatchField, volMesh> >&
    );


    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const DimensionedField<Type, volMesh>&
    );

    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const tmp<DimensionedField<Type, volMesh> >&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvcVolumeIntegrate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
