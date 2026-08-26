/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::LagrangianPatch

Description
    Base class for Lagrangian patches

SourceFiles
    LagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianPatch_H
#define LagrangianPatch_H

#include "polyPatch.H"
#include "LagrangianSubMesh.H"
#include "LagrangianFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianBoundaryMesh;
class objectRegistry;

/*---------------------------------------------------------------------------*\
                         Class LagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class LagrangianPatch
{
    // Private Data

        //- Reference to the corresponding polyPatch
        const polyPatch& poly_;

        //- Reference to boundary mesh
        const LagrangianBoundaryMesh& boundaryMesh_;

        //- The sub-mesh for this patch
        mutable autoPtr<LagrangianSubMesh> meshPtr_;


public:

    //- Runtime type information
    TypeName(polyPatch::typeName_());


    //- Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            LagrangianPatch,
            polyPatch,
            (
                const polyPatch& poly,
                const LagrangianBoundaryMesh& boundaryMesh
            ),
            (poly, boundaryMesh)
        );


    // Constructors

        //- Construct from a patch and a boundary mesh
        LagrangianPatch(const polyPatch&, const LagrangianBoundaryMesh&);

        //- Disallow default bitwise copy construction
        LagrangianPatch(const LagrangianPatch&) = delete;


    // Selectors

        //- Return a pointer to a new LagrangianPatch of specified type
        //  corresponding to the given polyPatch
        static autoPtr<LagrangianPatch> New
        (
            const word& type,
            const polyPatch& poly,
            const LagrangianBoundaryMesh& boundaryMesh
        );

        //- Return a pointer to a new LagrangianPatch corresponding to the
        //  given polyPatch
        static autoPtr<LagrangianPatch> New
        (
            const polyPatch& poly,
            const LagrangianBoundaryMesh& boundaryMesh
        );


    //- Destructor
    virtual ~LagrangianPatch();


    // Member Functions

        //- Return the polyPatch
        const polyPatch& poly() const
        {
            return poly_;
        }

        //- Return name
        virtual const word& name() const
        {
            return poly_.name();
        }

        //- Return whether or not this patch is coupled
        virtual bool coupled() const
        {
            return poly_.coupled();
        }

        //- Return the boundary mesh
        const LagrangianBoundaryMesh& boundaryMesh() const
        {
            return boundaryMesh_;
        }

        //- Return the sub-mesh associated with this patch
        virtual const LagrangianSubMesh& mesh() const;

        //- Return the local object registry
        const objectRegistry& db() const;

        //- Return time
        const Time& time() const;

        //- Lookup and return the patchField of the named field from the
        //  local objectRegistry.
        template<class GeometricField, class Type>
        const typename GeometricField::PatchField& lookupPatchField
        (
            const word& name
        );

        //- Initialise evaluation of changes in elements that have tracked to
        //  this patch
        virtual void initEvaluate
        (
            PstreamBuffers&,
            LagrangianMesh& mesh,
            const LagrangianInternalScalarDynamicField& fraction
        ) const;

        //- Evaluate changes in elements that have tracked to this patch
        virtual void evaluate
        (
            PstreamBuffers&,
            LagrangianMesh&,
            const LagrangianInternalScalarDynamicField& fraction
        ) const;

        //- Update for mesh changes
        virtual void partition() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianPatch&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
