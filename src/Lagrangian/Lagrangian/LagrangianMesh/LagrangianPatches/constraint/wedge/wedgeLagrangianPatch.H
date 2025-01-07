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
    Foam::wedgeLagrangianPatch

Description
    Wedge Lagrangian patch. This is used for the patches that have normals that
    point tangentially around an axisymmetric geometry.

SourceFiles
    wedgeLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeLagrangianPatch_H
#define wedgeLagrangianPatch_H

#include "LagrangianPatch.H"
#include "wedgePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class wedgeLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class wedgeLagrangianPatch
:
    public LagrangianPatch
{
private:

    // Private Data

        //- Reference to the wedge poly patch
        const wedgePolyPatch& wedgePatch_;

        //- Between crossFaces/evaluate and partition the elements have moved
        //  to the opposite patch, but the mesh has not been partitioned. So,
        //  the sub-mesh that relates to a wedge patch's elements is that of
        //  the opposite patch. This flag identifies when this is the case.
        mutable bool isOppositePatchMesh_;


public:

    //- Runtime type information
    TypeName(wedgePolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        wedgeLagrangianPatch(const polyPatch&, const LagrangianBoundaryMesh&);


    //- Destructor
    virtual ~wedgeLagrangianPatch();


    // Member Functions

        //- Return the wedge patch
        const wedgePolyPatch& wedgePatch() const
        {
            return wedgePatch_;
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
