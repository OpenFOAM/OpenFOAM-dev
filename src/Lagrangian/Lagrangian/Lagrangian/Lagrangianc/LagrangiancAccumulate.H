/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::Lagrangianc

Description
    Functions for accumulating Lagrangian values into cell fields

SourceFiles
    LagrangiancAccumulate.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangiancAccumulate_H
#define LagrangiancAccumulate_H

#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Namespace Lagrangianc functions Declaration
\*---------------------------------------------------------------------------*/

namespace Lagrangianc
{
    template<class CellMesh, class Type, template<class> class PrimitiveField>
    tmp<DimensionedField<Type, CellMesh>> accumulate
    (
        const DimensionedField<Type, LagrangianMesh, PrimitiveField>& lPsi,
        const word& name
    );

    template<class CellMesh, class Type, template<class> class PrimitiveField>
    tmp<DimensionedField<Type, CellMesh>> accumulate
    (
        const DimensionedField<Type, LagrangianMesh, PrimitiveField>& lPsi
    );

    template<class CellMesh, class Type, template<class> class PrimitiveField>
    void accumulate
    (
        const LagrangianSubField<Type, PrimitiveField>& lPsi,
        DimensionedField<Type, CellMesh>& vPsi
    );

    template<class CellMesh, class Type, template<class> class PrimitiveField>
    void accumulate
    (
        const LagrangianSubField<Type, PrimitiveField>& lPsi,
        DimensionedField<Type, CellMesh>& vPsi,
        const word& name
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangiancAccumulate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
