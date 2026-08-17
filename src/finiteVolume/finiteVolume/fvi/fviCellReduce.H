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
    Construct a volInternalField from a surface field using a combine operator

SourceFiles
    fviCellReduce.C

\*---------------------------------------------------------------------------*/

#ifndef fviCellReduce_H
#define fviCellReduce_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Namespace fvi functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvi
{
    template<class Type, class CombineOp>
    tmp<VolInternalField<Type>> cellReduce
    (
        const SurfaceField<Type>&,
        const CombineOp& cop,
        const Type& nullValue = pTraits<Type>::zero
    );

    template<class Type, class CombineOp>
    tmp<VolInternalField<Type>> cellReduce
    (
        const tmp<SurfaceField<Type>>&,
        const CombineOp& cop,
        const Type& nullValue = pTraits<Type>::zero
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fviCellReduce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
