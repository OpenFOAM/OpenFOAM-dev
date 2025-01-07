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
    Foam::cyclicLagrangianPatch

Description
    Cyclic Lagrangian patch. This is used for the patches that represent
    translated or rotated repetitions of the geometry.

SourceFiles
    cyclicLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicLagrangianPatch_H
#define cyclicLagrangianPatch_H

#include "LagrangianPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cyclicLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicLagrangianPatch
:
    public LagrangianPatch
{
private:

    // Private Data

        //- Reference to the cyclic poly patch
        const cyclicPolyPatch& cyclicPatch_;

        //- Between crossFaces/evaluate and partition the elements have moved
        //  to the neighbour patch, but the mesh has not been partitioned. So,
        //  the sub-mesh that relates to a cyclic patch's elements is that of
        //  the neighbour patch. This flag identifies when this is the case.
        mutable bool isNbrPatchMesh_;


public:

    //- Runtime type information
    TypeName(cyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        cyclicLagrangianPatch(const polyPatch&, const LagrangianBoundaryMesh&);


    //- Destructor
    virtual ~cyclicLagrangianPatch();


    // Member Functions

        //- Return the cyclic patch
        const cyclicPolyPatch& cyclicPatch() const
        {
            return cyclicPatch_;
        }

        //- Return a reference to the cyclic patch's transform
        const transformer& transform() const
        {
            return cyclicPatch_.transform();
        }

        //- Return the sub-mesh associated with this patch
        virtual const LagrangianSubMesh& mesh() const;

        //- Evaluate changes in elements that have tracked to this patch
        virtual void evaluate
        (
            PstreamBuffers&,
            LagrangianMesh&,
            const LagrangianScalarInternalDynamicField& fraction
        ) const;

        //- Update following partitioning of the mesh
        virtual void partition() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
