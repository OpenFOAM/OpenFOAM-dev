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

Class
    Foam::cellPointLagrangianAccumulator

Description
    Helper class for the cell-point Lagrangian accumulation scheme. This is a
    mesh object in order to provide a place to store re-usable accumulation
    weights and workspace.

SourceFiles
    cellPointLagrangianAccumulator.C

\*---------------------------------------------------------------------------*/

#ifndef cellPointLagrangianAccumulator_H
#define cellPointLagrangianAccumulator_H

#include "CompactListList.H"
#include "DemandDrivenMeshObject.H"
#include "LagrangianSubFieldsFwd.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cellPointLagrangianAccumulator Declaration
\*---------------------------------------------------------------------------*/

class cellPointLagrangianAccumulator
:
    public DemandDrivenMeshObject
    <
        polyMesh,
        MoveableMeshObject,
        cellPointLagrangianAccumulator
    >
{
    // Private Data

        //- Weights used to accumulate the point-field back into the cell-field
        CompactListList<scalar> pointCellWeights_;

        //- Map from point to accumulating point
        mutable labelList pointAccumulatingPoint_;

        //- Map from accumulating point to point
        mutable DynamicList<label> accumulatingPointPoint_;

        //- Accumulated values on each accumulating point
        #define DECLARE_ACCUMULATING_POINT_TYPES(Type, nullArg)       \
            mutable autoPtr<DynamicList<Type>>                        \
                CAT3(accumulatingPoint, CAPITALIZE(Type), ValuesPtr_);
        FOR_ALL_FIELD_TYPES(DECLARE_ACCUMULATING_POINT_TYPES);
        #undef DECLARE_ACCUMULATING_POINT_TYPES


    // Private Member Functions

        //- Calculate the point-cell weights
        void calcPointCellWeights();

        //- Access accumulated values
        template<class Type>
        DynamicList<Type>& accumulatingPointValues() const;


protected:

    // Protected Constructors

        //- Allow the base mesh object to use the protected constructor
        friend class DemandDrivenMeshObject
        <
            polyMesh,
            MoveableMeshObject,
            cellPointLagrangianAccumulator
        >;

        //- Construct from mesh
        explicit cellPointLagrangianAccumulator(const polyMesh& mesh);


public:

    // Static Data Members

        //- Runtime type information
        TypeName("cellPointLagrangianAccumulator");


    // Constructors

        //- Disallow default bitwise copy construction
        cellPointLagrangianAccumulator(const cellPointLagrangianAccumulator&) =
            delete;


    //- Destructor
    ~cellPointLagrangianAccumulator();


    // Member Functions

        //- Update for mesh motion
        virtual bool movePoints();

        //- Accumulate the Lagrangian field into a primitive cell field
        template<class Type>
        void accumulate
        (
            const LagrangianSubSubField<Type>& lPsi,
            Field<Type>& cPsi
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cellPointLagrangianAccumulator&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cellPointLagrangianAccumulatorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
