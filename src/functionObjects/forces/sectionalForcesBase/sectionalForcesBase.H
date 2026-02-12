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

Class
    Foam::functionObjects::sectionalForcesBase

Description
    Base class for sectional forces function objects

SourceFiles
    sectionalForcesBase.C

\*---------------------------------------------------------------------------*/

#ifndef sectionalForcesBase_functionObject_H
#define sectionalForcesBase_functionObject_H

#include "fvMeshFunctionObject.H"
#include "indirectPrimitivePatch.H"
#include "patchCutPlot.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class sectionalForcesBase Declaration
\*---------------------------------------------------------------------------*/

class sectionalForcesBase
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Mesh patches to integrate the forces over
        labelHashSet patchSet_;

        //- Primitive patch containing all mesh patches
        mutable autoPtr<indirectPrimitivePatch> patchPtr_;

        //- Name of pressure field
        word pName_;

        //- Name of velocity field
        word UName_;

        //- Name of density field (optional)
        word rhoName_;

        //- The name of the phase (optional)
        word phaseName_;

        //- Reference density needed for incompressible calculations
        scalar rhoRef_;

        //- Reference pressure
        scalar pRef_;

        //- Weights
        mutable autoPtr<List<patchCutPlot::weight>> weightsPtr_;


    // Private Member Functions

        //- Multiply the argument by the volume fraction (if any)
        template<class Type>
        tmp<SurfaceField<Type>> timesAlpha
        (
            const tmp<SurfaceField<Type>> psi
        ) const;

        //- Multiply the argument by the volume fraction (if any)
        template<class Type>
        tmp<Field<Type>> timesAlpha
        (
            const tmp<Field<Type>> psi,
            const label patchi
        ) const;

        //- Multiply the argument by the density
        template<class Type>
        tmp<SurfaceField<Type>> timesRho
        (
            const tmp<SurfaceField<Type>> psi
        ) const;

        //- Multiply the argument by the volume fraction (if any) and density
        template<class Type>
        tmp<SurfaceField<Type>> timesAlphaRho
        (
            const tmp<SurfaceField<Type>> psi
        ) const;

        //- Return the (dynamic) pressure
        tmp<volScalarField> p() const;

        //- Return the surface stress tensor
        tmp<surfaceVectorField> devTau() const;


protected:

    // Protected Member Functions

        //- Clear the cached weights
        virtual void clear();

        //- Access the primitive patch
        const indirectPrimitivePatch& patch() const;

        //- Clear the patch
        void clearPatch();

        //- Clear the patch geometry
        void clearPatchGeom();

        //- Return the normal to the section planes
        virtual vector normal() const = 0;

        //- Return the origin of the plot. This can be an arbitrary point
        //  somewhere along the 1-D 'beam' that the geometry is reduced to.
        virtual point origin() const = 0;

        //- Return the distances from the origin to the section planes. This is
        //  the coordinate that we plot against.
        virtual tmp<scalarField> distances() const = 0;

        //- Return the distance from the origin to the patch points
        //  in the normal direction
        tmp<scalarField> patchPointDistances() const;

        //- Access the weights
        const List<patchCutPlot::weight>& weights() const;

        //- Return the output path
        fileName outputPath() const;

        //- Calculate the forces and moments at the cuts
        void addFluid
        (
            vectorField& force,
            vectorField& moment
        ) const;


public:

    //- Runtime type information
    TypeName("sectionalForcesBase");


    // Constructors

        //- Construct from Time and dictionary
        sectionalForcesBase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sectionalForcesBase(const sectionalForcesBase&) = delete;


    //- Destructor
    virtual ~sectionalForcesBase();


    // Member Functions

        //- Read the sectionalForcesBase data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh&);

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sectionalForcesBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
