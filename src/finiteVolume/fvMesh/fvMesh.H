/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::fvMesh

Description
    Mesh data needed to do the Finite Volume discretisation.

    NOTE ON USAGE:
    fvMesh contains all the topological and geometric information
    related to the mesh.  It is also responsible for keeping the data
    up-to-date.  This is done by deleting the cell volume, face area,
    cell/face centre, addressing and other derived information as
    required and recalculating it as necessary.  The fvMesh therefore
    reserves the right to delete the derived information upon every
    topological (mesh refinement/morphing) or geometric change (mesh
    motion).  It is therefore unsafe to keep local references to the
    derived data outside of the time loop.

SourceFiles
    fvMesh.C
    fvMeshGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef fvMesh_H
#define fvMesh_H

#include "polyMesh.H"
#include "lduMesh.H"
#include "primitiveMesh.H"
#include "fvBoundaryMesh.H"
#include "surfaceInterpolation.H"
#include "fvSchemes.H"
#include "fvSolution.H"
#include "data.H"
#include "DimensionedField.H"
#include "SlicedDimensionedField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "slicedVolFieldsFwd.H"
#include "slicedSurfaceFieldsFwd.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshLduAddressing;
class fvMeshTopoChanger;
class fvMeshDistributor;
class fvMeshMover;
class volMesh;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class fvMesh Declaration
\*---------------------------------------------------------------------------*/

class fvMesh
:
    public polyMesh,
    public lduMesh,
    public surfaceInterpolation,
    public fvSchemes,
    public fvSolution,
    public data
{
    // Private Data

        //- Boundary mesh
        fvBoundaryMesh boundary_;

        //- The topo-changer function class
        autoPtr<fvMeshTopoChanger> topoChanger_;

        //- The distributor function class
        autoPtr<fvMeshDistributor> distributor_;

        //- The mover function class
        autoPtr<fvMeshMover> mover_;


    // Demand-driven data

        mutable fvMeshLduAddressing* lduPtr_;

        //- Current time index for cell volumes
        //  Note.  The whole mechanism will be replaced once the
        //  dimensionedField is created and the dimensionedField
        //  will take care of the old-time levels.
        mutable label curTimeIndex_;

        //- Cell volumes old time level
        mutable SlicedDimensionedField<scalar, volMesh>* VPtr_;

        //- Cell volumes old time level
        mutable DimensionedField<scalar, volMesh>* V0Ptr_;

        //- Cell volumes old-old time level
        mutable DimensionedField<scalar, volMesh>* V00Ptr_;

        //- Face area vectors
        mutable slicedSurfaceVectorField* SfPtr_;

        //- Mag face area vectors
        mutable slicedSurfaceScalarField* magSfPtr_;

        //- Cell centres
        mutable slicedVolVectorField* CPtr_;

        //- Face centres
        mutable slicedSurfaceVectorField* CfPtr_;

        //- Face motion fluxes
        mutable surfaceScalarField* phiPtr_;


    // Private Member Functions

        // Storage management

            //- Clear geometry but not the old-time cell volumes
            void clearGeomNotOldVol();

            //- Clear geometry like clearGeomNotOldVol but recreate any
            //  geometric demand-driven data that was set
            void updateGeomNotOldVol();

            //- Clear geometry
            void clearGeom();

            //- Clear addressing
            void clearAddressing(const bool isMeshUpdate = false);

            //- Preserve old volume(s)
            void storeOldVol(const scalarField&);


       // Make geometric data

            void makeSf() const;
            void makeMagSf() const;

            void makeC() const;
            void makeCf() const;


public:

    // Public Typedefs

        typedef fvMesh Mesh;
        typedef fvBoundaryMesh BoundaryMesh;


    // Declare name of the class and its debug switch
    ClassName("fvMesh");


    // Constructors

        //- Construct from IOobject
        //  with the option to not instantiate the mesh changers
        explicit fvMesh(const IOobject& io, const bool changers = true);

        //- Construct from cellShapes with boundary.
        fvMesh
        (
            const IOobject& io,
            pointField&& points,
            const cellShapeList& shapes,
            const faceListList& boundaryFaces,
            const wordList& boundaryPatchNames,
            const PtrList<dictionary>& boundaryDicts,
            const word& defaultBoundaryPatchName,
            const word& defaultBoundaryPatchType,
            const bool syncPar = true
        );

        //- Construct from components without boundary.
        //  Boundary is added using addFvPatches() member function
        fvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            labelList&& allOwner,
            labelList&& allNeighbour,
            const bool syncPar = true
        );

        //- Construct without boundary from cells rather than owner/neighbour.
        //  Boundary is added using addPatches() member function
        fvMesh
        (
            const IOobject& io,
            pointField&& points,
            faceList&& faces,
            cellList&& cells,
            const bool syncPar = true
        );

        //- Disallow default bitwise copy construction
        fvMesh(const fvMesh&);


    //- Destructor
    virtual ~fvMesh();


    // Member Functions

        // Helpers

            //- Add boundary patches. Constructor helper
            void addFvPatches
            (
                const List<polyPatch*>&,
                const bool validBoundary = true
            );

            //- Update the mesh based on the mesh files saved in time
            //  directories
            virtual readUpdateState readUpdate();


        // Access

            //- Return the top-level database
            const Time& time() const
            {
                return polyMesh::time();
            }

            //- Return the object registry - resolve conflict polyMesh/lduMesh
            virtual const objectRegistry& thisDb() const
            {
                return polyMesh::thisDb();
            }

            //- Return reference to name
            //  Note: name() is currently ambiguous due to derivation from
            //  surfaceInterpolation
            const word& name() const
            {
                return polyMesh::name();
            }

            //- Return reference to boundary mesh
            const fvBoundaryMesh& boundary() const;

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return boundary().interfaces();
            }

            //- Return communicator used for parallel communication
            virtual label comm() const
            {
                return polyMesh::comm();
            }

            //- Internal face owner
            const labelUList& owner() const
            {
                return lduAddr().lowerAddr();
            }

            //- Internal face neighbour
            const labelUList& neighbour() const
            {
                return lduAddr().upperAddr();
            }

            //- Return the topo-changer function class
            const fvMeshTopoChanger& topoChanger() const;

            //- Return the distributor function class
            const fvMeshDistributor& distributor() const;

            //- Return the mover function class
            const fvMeshMover& mover() const;

            //- Return cell volumes
            const DimensionedField<scalar, volMesh>& V() const;

            //- Return old-time cell volumes
            const DimensionedField<scalar, volMesh>& V0() const;

            //- Access old-time cell volumes
            DimensionedField<scalar, volMesh>& V0Ref();

            //- Return old-old-time cell volumes
            const DimensionedField<scalar, volMesh>& V00() const;

            //- Access old-old-time cell volumes
            DimensionedField<scalar, volMesh>& V00Ref();

            //- Return sub-cycle cell volumes
            tmp<DimensionedField<scalar, volMesh>> Vsc() const;

            //- Return sub-cycle old-time cell volumes
            tmp<DimensionedField<scalar, volMesh>> Vsc0() const;

            //- Return cell face area vectors
            const surfaceVectorField& Sf() const;

            //- Return cell face area magnitudes
            const surfaceScalarField& magSf() const;

            //- Return cell centres as volVectorField
            const volVectorField& C() const;

            //- Return face centres as surfaceVectorField
            const surfaceVectorField& Cf() const;

            //- Return face deltas as surfaceVectorField
            tmp<surfaceVectorField> delta() const;

            //- Return cell face motion fluxes
            const surfaceScalarField& phi() const;

            //- Access cell face motion fluxes
            surfaceScalarField& phiRef();

            //- Return a labelType of valid component indicators
            //  1 : valid (solved)
            // -1 : invalid (not solved)
            template<class Type>
            typename pTraits<Type>::labelType validComponents() const;


        // Edit

            //- Is mesh dynamic
            virtual bool dynamic() const;

            //- Update the mesh for both mesh motion and topology change
            virtual bool update();

            //- Clear all geometry and addressing
            void clearOut();

            //- Update mesh corresponding to the given map
            virtual void updateMesh(const mapPolyMesh& map);

            //- Update mesh corresponding to the given distribution map
            virtual void distribute(const mapDistributePolyMesh& map);

            //- Move points, returns volumes swept by faces in motion
            virtual tmp<scalarField> movePoints(const pointField&);

            //- Map all fields in time using given map.
            virtual void mapFields(const mapPolyMesh& mpm);

            //- Add/insert single patch. If validBoundary the new situation
            //  is consistent across processors.
            virtual void addPatch
            (
                const label insertPatchi,
                const polyPatch& patch,
                const dictionary& patchFieldDict,
                const word& defaultPatchFieldType,
                const bool validBoundary
            );

            //- Reorder and trim existing patches. If validBoundary the new
            //  situation is consistent across processors
            virtual void reorderPatches
            (
                const labelUList& newToOld,
                const bool validBoundary
            );

            //- Remove boundary patches. Warning: fvPatchFields hold ref to
            //  these fvPatches.
            void removeFvBoundary();


        // Write

            //- Write the underlying polyMesh and other data
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write = true
            ) const;

            //- Write mesh using IO settings from time
            virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMesh&);

        bool operator!=(const fvMesh&) const;
        bool operator==(const fvMesh&) const;
};


template<>
typename pTraits<sphericalTensor>::labelType
fvMesh::validComponents<sphericalTensor>() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
