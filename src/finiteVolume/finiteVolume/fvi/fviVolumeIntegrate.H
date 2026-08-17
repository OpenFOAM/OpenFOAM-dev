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
    Volume integrate volField creating a volField.

    Volume integrate volField over the whole domain
    creating a dimensioned\<Type\>

SourceFiles
    fviVolumeIntegrate.C

\*---------------------------------------------------------------------------*/


#ifndef fviVolumeIntegrate_H
#define fviVolumeIntegrate_H

#include "volFieldsFwd.H"
#include "primitiveFieldsFwd.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvi functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvi
{
    template<class Type>
    tmp<Field<Type>> volumeIntegrate(const VolField<Type>&);

    template<class Type>
    tmp<Field<Type>> volumeIntegrate(const tmp<VolField<Type>>&);


    template<class Type>
    tmp<Field<Type>> volumeIntegrate(const DimensionedField<Type, fvMesh>&);

    template<class Type>
    tmp<Field<Type>> volumeIntegrate
    (
        const tmp<DimensionedField<Type, fvMesh>>&
    );


    template<class Type>
    dimensioned<Type> domainIntegrate(const VolField<Type>&);

    template<class Type>
    dimensioned<Type> domainIntegrate(const tmp<VolField<Type>>&);


    template<class Type>
    dimensioned<Type> domainIntegrate(const DimensionedField<Type, fvMesh>&);

    template<class Type>
    dimensioned<Type> domainIntegrate
    (
        const tmp<DimensionedField<Type, fvMesh>>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fviVolumeIntegrate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
