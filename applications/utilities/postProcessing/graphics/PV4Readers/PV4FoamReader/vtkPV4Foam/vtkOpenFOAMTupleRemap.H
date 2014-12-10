/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    vtkPV4Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkOpenFOAMTupleRemap_H
#define vtkOpenFOAMTupleRemap_H

// OpenFOAM includes
#include "StaticAssert.H"
#include "Swap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void vtkOpenFOAMTupleRemap(float vec[]);


// a symmTensor specialization to remap OpenFOAM -> ParaView naming order
// Qt/Core/pqScalarBarRepresentation.cxx defines this order
//     { "XX", "YY", "ZZ", "XY", "YZ", "XZ" }
// in pqScalarBarRepresentation::getDefaultComponentLabel()
// whereas OpenFOAM uses this order
//     { XX, XY, XZ, YY, YZ, ZZ }
//
// for extra safety, assert that symmTensor indeed has 6 components
StaticAssert(Foam::symmTensor::nComponents == 6);


// Template specialization for symmTensor
template<>
inline void vtkOpenFOAMTupleRemap<Foam::symmTensor>(float vec[])
{
    Foam::Swap(vec[1], vec[3]);   // swap XY <-> YY
    Foam::Swap(vec[2], vec[5]);   // swap XZ <-> ZZ
}


template<class Type>
inline void vtkOpenFOAMTupleRemap(float vec[])
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
