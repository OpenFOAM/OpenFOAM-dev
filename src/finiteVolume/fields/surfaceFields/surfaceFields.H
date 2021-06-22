/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

InClass
    Foam::surfaceFields

Description
    Foam::surfaceFields

SourceFiles
    surfaceFields.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFields_H
#define surfaceFields_H

#include "GeometricFields.H"
#include "surfaceMesh.H"
#include "fvMesh.H"
#include "fvsPatchFields.H"
#include "surfaceFieldsFwd.H"
#include "calculatedFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Check if surfaceField is a flux
//  i.e. an extensive field oriented with the face requiring
//  flipping after mapping.
//  Currently it is assumed that only surfaceScalarFields can be fluxes
template<class Type>
inline bool isFlux(const DimensionedField<Type, surfaceMesh>& df)
{
    return false;
}

//- Check if surfaceScalarField is a flux
//  i.e. an extensive field oriented with the face requiring
//  flipping after mapping.
//  Currently it is assumed that fluxes are either volumetric or mass
template<>
inline bool isFlux<scalar>(const DimensionedField<scalar, surfaceMesh>& df)
{
    if (df.dimensions() == dimFlux || df.dimensions() == dimMassFlux)
    {
        return true;
    }
    else
    {
        return false;
    }
}

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
