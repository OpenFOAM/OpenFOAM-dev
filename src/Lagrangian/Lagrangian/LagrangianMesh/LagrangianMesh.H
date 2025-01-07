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
    Foam::LagrangianMesh

Description
    Class containing Lagrangian geometry and topology

SourceFiles
    LagrangianMesh.C
    LagrangianMeshTemplates.C
    LagrangianMeshI.H
    LagrangianMeshLocation.H

\*---------------------------------------------------------------------------*/

#ifndef LagrangianMesh_H
#define LagrangianMesh_H

#include "GeoMesh.H"
#include "barycentricField.H"
#include "barycentricIODynamicField.H"
#include "LagrangianBoundaryMesh.H"
#include "LagrangianFieldsFwd.H"
#include "LagrangianSchemes.H"
#include "LagrangianSolution.H"
#include "LagrangianState.H"
#include "LagrangianSubMesh.H"
#include "labelIODynamicField.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianInjection;
class nonConformalCyclicPolyPatch;

/*---------------------------------------------------------------------------*\
                       Class LagrangianMesh Declaration
\*---------------------------------------------------------------------------*/

class LagrangianMesh
:
    public objectRegistry,
    public GeoMesh<polyMesh>
{
public:

    // Public Enumerations

        //- Enumeration of the permutation algorithm
        enum class permutationAlgorithm
        {
            copy,
            inPlace
        };

        //- Permutation algorithm names
        static const NamedEnum<permutationAlgorithm, 2>
            permutationAlgorithmNames_;

        //- Enumeration of the partitioning algorithm
        enum class partitioningAlgorithm
        {
            bin,
            quick
        };

        //- Partitioning algorithm names
        static const NamedEnum<partitioningAlgorithm, 2>
            partitioningAlgorithmNames_;

        //- Enumeration for the locations of searched positions
        enum class location
        {
            inCell,
            onBoundary,
            outsideMesh
        };


private:

    // Private Data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Local coordinates of the elements within their tetrahedra
        barycentricIODynamicField coordinates_;

        //- Indices of the cells that contain the tetrahedra
        labelIODynamicField celli_;

        //- Indices of the faces that form the bases of the tetrahedra
        labelIODynamicField facei_;

        //- Indices of the face-triangles that form the bases of the tetrahedra
        labelIODynamicField faceTrii_;

        //- Boundary mesh
        LagrangianBoundaryMesh boundary_;

        //- A special sub-mesh which spans the entire mesh
        LagrangianSubMesh subAll_;

        //- The elements' states
        autoPtr<DynamicList<LagrangianState>> statesPtr_;

        //- Offsets defining ranges associated with the elements' groups
        autoPtr<labelList> offsetsPtr_;

        //- For every patch, a list of associated non-conformal cyclic patch
        //  indices. Empty if there are no non-conformal cyclic patches.
        autoPtr<labelListList> origPatchNccPatchisPtr_;

        //- For every patch, a list of associated non-conformal cyclic patch
        //  pointers. Empty if there are no non-conformal cyclic patches.
        autoPtr<List<UPtrList<const nonConformalCyclicPolyPatch>>>
            origPatchNccPatchesPtr_;

        //- For every non-conformal cyclic patch, the indices of the associated
        //  non-conformal processor-cyclic patches. Empty if there are no
        //  non-conformal cyclic patches.
        autoPtr<labelListList> nccPatchProcNccPatchisPtr_;

        //- The receiving patch face for non-conformal cyclic transfers
        autoPtr<DynamicList<label>> receivePatchFacePtr_;

        //- The receiving position for non-conformal cyclic transfers
        autoPtr<DynamicList<point>> receivePositionPtr_;

        //- The step fraction less than the maximum reached so far. See
        //  tracking.H for details.
        autoPtr<LagrangianDynamicField<scalar>> fractionBehindPtr_;

        //- The number of tracks carried out that ended in a step fraction less
        //  than the maximum reached so far. See tracking.H for details.
        autoPtr<LagrangianDynamicField<label>> nTracksBehindPtr_;

        //- Sub-mesh index
        mutable label subMeshIndex_;

        //- Schemes created on demand
        mutable autoPtr<LagrangianSchemes> schemesPtr_;

        //- Solution controls created on demand
        mutable autoPtr<LagrangianSolution> solutionPtr_;


    // Private Member Functions

        // Checks

            //- Check that the field is the same size as the mesh
            template<class Field>
            inline void checkFieldSize(const Field& field) const;

            //- Check that a pointer is valid
            template<class Type>
            inline void checkPtr
            (
                const autoPtr<Type>& ptr,
                const word& name
            ) const;

            //- Print a row of the group table
            void printGroups(const bool header) const;


        // Partitioning and Permutation

            //- Partition the states using the bin-sort-based algorithm
            labelList partitionBin
            (
                labelList& offsets,
                const List<LagrangianState>& states
            ) const;

            //- Partition the states using a quick-sort-based algorithm
            labelList partitionQuick
            (
                labelList& offsets,
                const List<LagrangianState>& states
            ) const;

            //- Reorder and resize all registered fields using the given
            //  permutation
            void permuteAndResizeFields(const labelList& permutation);

            //- Reorder a list with the given permutation
            template<class Type>
            static void permuteList
            (
                const labelList& permutation,
                UList<Type>& list
            );

            //- Reorder a list using the copy algorithm
            template<class Type>
            static void permuteListCopy
            (
                const labelList& permutation,
                UList<Type>& list
            );

            //- Reorder a list using the in-place algorithm
            template<class Type>
            static void permuteListInPlace
            (
                const labelList& permutation,
                UList<Type>& list
            );

            //- Resize a container to match the mesh
            template<class Container>
            void resizeContainer(Container& container) const;


        // Addition

            //- Append specified elements in the mesh with the given geometry
            //  and topology. Does nothing to fields.
            LagrangianSubMesh append
            (
                const barycentricField& coordinates,
                const labelField& celli,
                const labelField& facei,
                const labelField& faceTrii
            );

            //- Birth specified elements in the mesh from the given list of
            //  parent indices. Does nothing to fields.
            LagrangianSubMesh append(const labelList& parents);

            //- Expand and set values for explicitly specified fields following
            //  an injection or birth event
            template<class Type, template<class> class GeoField>
            void appendSpecifiedField
            (
                const LagrangianSubMesh& appendMesh,
                GeoField<Type>& geoField,
                const Field<Type>& field
            ) const;

            //- Expand and set values for explicitly specified fields following
            //  an injection or birth event
            template<class Type, template<class> class GeoField>
            bool appendSpecifiedField
            (
                const LagrangianSubMesh& appendMesh,
                const word& fieldName,
                const Field<Type>& field
            ) const;

            //- Expand and set values for explicitly specified fields following
            //  an injection or birth event
            template<class Type, class ... FieldNamesAndFields>
            wordHashSet appendSpecifiedFields
            (
                const LagrangianSubMesh& appendMesh,
                const word& fieldName,
                const Field<Type>& field,
                const FieldNamesAndFields& ... fieldNamesAndFields
            ) const;

            //- Termination clause for the above
            wordHashSet appendSpecifiedFields
            (
                const LagrangianSubMesh& appendMesh
            ) const;

            //- Expand and set values for unspecified fields following an
            //  injection
            void injectUnspecifiedFields
            (
                const LagrangianSubMesh& injectionMesh,
                const wordHashSet& specifiedFieldNames
            );

            //- Expand and set values for unspecified fields following an
            //  injection
            void injectUnspecifiedFields
            (
                const LagrangianInjection& injection,
                const LagrangianSubMesh& injectionMesh,
                const wordHashSet& specifiedFieldNames
            );

            //- Expand and set values for unspecified fields following a birth
            void birthUnspecifiedFields
            (
                const labelList& parents,
                const LagrangianSubMesh& birthMesh,
                const wordHashSet& specifiedFieldNames
            );


public:

    // Friend classes

        //- The non-conformal cyclic patch needs to access the receiving
        //  information populated during tracking
        friend class nonConformalCyclicLagrangianPatch;

        //- The processor patch communicates and adds elements to the receiving
        //  side, but field values are added later by the corresponding patch
        //  field. This is not permitted by the public interface in which
        //  elements and associated field values are required to be added
        //  simultaneously so that they are guaranteed to remain consistent.
        //  The requirements of the processor patch and patch fields are
        //  considered a special case, so access to the private append methods
        //  is permitted to these classes.
        friend class processorLagrangianPatch;

        //- See above
        template<class Type>
        friend class processorLagrangianPatchField;

        //- See above
        friend class nonConformalProcessorCyclicLagrangianPatch;

        //- See above
        template<class Type>
        friend class nonConformalProcessorCyclicLagrangianPatchField;


    // Public classes

        //- Class to define the scope of Lagrangian mesh state changes
        class changer
        {
            // Private Data

                //- Reference to the Lagrangian mesh
                LagrangianMesh& mesh_;


            // Private Member Functions

                //- Construct the mesh's non-conformal data
                void constructNonConformal() const;

                //- Construct the mesh's behind data
                void constructBehind() const;


        public:

            // Constructors

                //- Construct for a Lagrangian mesh with a given state
                changer
                (
                    LagrangianMesh& mesh,
                    const LagrangianState state
                );

                //- Construct for a Lagrangian mesh with given states
                changer
                (
                    LagrangianMesh& mesh,
                    const List<LagrangianState>& state
                );


            // Destructor
            ~changer();
        };

        //- Class to hold and index into the field references associated with a
        //  linear displacement
        class linearDisplacement
        {
            // Private Member Data

                //- Reference to the linear displacement
                const LagrangianSubVectorField& linear_;

        public:

            // Constructors

                //- Construct from a reference to the displacement
                linearDisplacement(const LagrangianSubVectorField& linear);


            // Member Operators

                //- The displacement for a given index
                inline const vector& operator()(const label i) const;

                //- The displacement for a given index and remaining fraction
                inline const vector& operator()
                (
                    const label i,
                    const scalar f
                ) const;
        };

        //- Class to hold and index into the field references associated with a
        //  parabolic displacement
        class parabolicDisplacement
        {
            // Private Member Data

                //- Reference to the linear displacement
                const LagrangianSubVectorField& linear_;

                //- Reference to the quadratic displacement
                const LagrangianSubVectorField& quadratic_;


        public:

            // Constructors

                //- Construct from references to the displacements
                parabolicDisplacement
                (
                    const LagrangianSubVectorField& linear,
                    const LagrangianSubVectorField& quadratic
                );


            // Member Operators

                //- The displacements for a given index
                inline Pair<vector> operator()(const label i) const;

                //- The displacement for a given index and remaining fraction
                inline vector operator()
                (
                    const label i,
                    const scalar f
                ) const;
        };

        //- Class to hold element-group indices, and associate the group
        //  indices with a given enumeration. Binary compatible with labelPair
        //  so that a UList<elementGroup<Enumeration>> can be cast to a
        //  UList<labelPair>. This allows the labelPair-based partition
        //  function to be used with an arbitrary enumeration without
        //  templating the implementation.
        template<class Enumeration>
        class elementGroup
        :
            private labelPair
        {
            //- Reference class returned by the group accessors
            class EnumerationRef
            {
            private:

                label& i_;

            public:

                EnumerationRef(label& i)
                :
                    i_(i)
                {}

                operator Enumeration() const
                {
                    return static_cast<Enumeration>(i_);
                }

                void operator=(const Enumeration e)
                {
                    i_ = static_cast<label>(e);
                }
            };

            //- Modify the element index
            label& element()
            {
                return first();
            }

            //- Access the element index
            label element() const
            {
                return first();
            }

            //- Modify the group enumeration
            EnumerationRef group()
            {
                return EnumerationRef(second());
            }

            //- Access the group enumeration
            Enumeration group() const
            {
                return static_cast<Enumeration>(second());
            }
        };


    // Public Static Data

        //- Run-time type information
        ClassName("LagrangianMesh");

        //- Instance prefix
        static const word prefix;

        //- Name of the coordinates field
        static const word coordinatesName;

        //- Name of the position field
        static const word positionName;

        //- Name of the state field
        static const word stateName;

        //- Name of the tracked fraction field
        static const word fractionName;

        //- Permutation algorithm
        static permutationAlgorithm permutationAlgorithm_;

        //- Partitioning algorithm
        static partitioningAlgorithm partitioningAlgorithm_;


    // Public Type Definitions

        //- Mesh type
        typedef LagrangianMesh Mesh;

        //- Boundary mesh type
        typedef LagrangianBoundaryMesh BoundaryMesh;

        //- Patch field type
        template<class Type>
        using PatchField = LagrangianPatchField<Type>;

        //- Field source type
        template<class Type>
        using FieldSource = LagrangianFieldSource<Type>;


    // Constructors

        //- Construct from a mesh and a name
        LagrangianMesh
        (
            const polyMesh& mesh,
            const word& name,
            const IOobject::readOption readOption=IOobject::READ_IF_PRESENT,
            const IOobject::writeOption writeOption=IOobject::AUTO_WRITE
        );

        //- Construct from a mesh and a name and a list of patch types
        LagrangianMesh
        (
            const polyMesh& mesh,
            const word& name,
            const wordList& wantedPatchTypes,
            const IOobject::readOption readOption=IOobject::READ_IF_PRESENT,
            const IOobject::writeOption writeOption=IOobject::AUTO_WRITE
        );

        //- Disallow default bitwise copy construction
        LagrangianMesh(const LagrangianMesh&) = delete;


    //- Destructor
    ~LagrangianMesh();


    // Member Functions

        // Access

            using objectRegistry::time;
            using objectRegistry::thisDb;
            using objectRegistry::operator!=;
            using objectRegistry::operator==;

            //- Access the mesh
            inline const polyMesh& mesh() const;

            //- Access the coordinates
            inline const barycentricIODynamicField& coordinates() const;

            //- Access the cell indices
            inline const labelIODynamicField& celli() const;

            //- Access the cell-face indices
            inline const labelIODynamicField& facei() const;

            //- Access the face-tet indices
            inline const labelIODynamicField& faceTrii() const;

            //- Return reference to boundary mesh
            inline const LagrangianBoundaryMesh& boundary() const;

            //- Return whether or not the mesh is changing
            inline bool changing() const;

            //- Return the states
            inline const List<LagrangianState>& states() const;

            //- Access the states
            inline List<LagrangianState>& states();

            //- Return the state for an element of the mesh, or a none state
            //  if not changing
            inline LagrangianState state(const label i) const;

            //- Return the state for an element of a sub mesh, or a none state
            //  if not changing
            inline LagrangianState state
            (
                const LagrangianSubMesh& subMesh,
                const label subi
            ) const;

            //- Get a sub-mesh index
            inline label subMeshIndex() const;

            //- Access the schemes
            const LagrangianSchemes& schemes() const;

            //- Access the solution controls
            const LagrangianSolution& solution() const;

            //- Lookup all current-time fields of the given type
            template<class GeoField>
            HashTable<const GeoField*> lookupCurrentFields
            (
                const bool strict = false
            ) const;

            //- Lookup all current-time fields of the given type
            template<class GeoField>
            HashTable<GeoField*> lookupCurrentFields
            (
                const bool strict = false
            );


        // Check

            //- Get the number of elements
            inline label size() const;

            //- Get the global number of elements
            inline label globalSize() const;

            //- Return size
            static inline label size(const LagrangianMesh& mesh);

            //- Return the number of states
            inline label nStates() const;

            //- Return the number of groups
            inline label nGroups() const;

            //- Convert a state to a group index
            inline label stateToGroupi(const LagrangianState state) const;

            //- Return a sub-mesh for the given group
            inline LagrangianSubMesh sub(const LagrangianGroup group) const;

            //- Return a sub-mesh for no elements
            inline LagrangianSubMesh subNone() const;

            //- Return a sub-mesh for all elements
            inline const LagrangianSubMesh& subAll() const;

            //- Return a sub-mesh for all incomplete elements
            inline LagrangianSubMesh subIncomplete() const;

            //- Return the global sizes of all the sub-meshes. A value of -1
            //  will be set for processor-patch sub-meshes.
            labelList subMeshGlobalSizes() const;

            //- Return the global positions
            tmp<LagrangianVectorInternalField> position() const;

            //- Return the global position of an element
            point position(const label i) const;

            //- Convert a position into a set of coordinates and a
            //  corresponding tetrahedron. Return a status flag indicating
            //  where the point is relative to the bounds of the mesh.
            location locate
            (
                const point& position,
                barycentric& coordinates,
                label& celli,
                label& facei,
                label& faceTrii,
                const scalar fraction
            ) const;

            //- Convert set of positions into a set of coordinates and a
            //  corresponding tetrahedron. Return status flags indicating
            //  where the points are relative to the bounds of the mesh.
            List<location> locate
            (
                const List<point>& position,
                List<barycentric>& coordinates,
                labelList& celli,
                labelList& facei,
                labelList& faceTrii,
                const scalarList& fraction
            ) const;


        // Modify

            //- Partition the mesh such that the groups are contiguous in memory
            void partition();

            //- Track the positions along the given displacements
            template<class Displacement>
            void track
            (
                const List<LagrangianState>& endState,
                const Displacement& displacement,
                const LagrangianSubScalarField& deltaFraction,
                LagrangianSubScalarSubField& fraction
            );

            //- Cross the faces
            void crossFaces
            (
                const LagrangianScalarInternalDynamicField& fraction
            );

            //- Inject specified elements into the mesh. This method does not
            //  know how to set values for the fields. If any fields are
            //  registered then this function will error.
            template<class ... FieldNamesAndFields>
            LagrangianSubMesh inject
            (
                const barycentricField& coordinates,
                const labelField& celli,
                const labelField& facei,
                const labelField& faceTrii,
                const FieldNamesAndFields& ... fieldNamesAndFields
            );

            //- Inject specified elements into the mesh. Fields are expanded
            //  accordingly. New field values can be specified manually by
            //  providing alternating name and field arguments. Source
            //  conditions will be used to set values for all other fields.
            //  Internal fields without source conditions must be specified
            //  manually or an error will result.
            template<class ... FieldNamesAndFields>
            LagrangianSubMesh inject
            (
                const LagrangianInjection& injection,
                const barycentricField& coordinates,
                const labelField& celli,
                const labelField& facei,
                const labelField& faceTrii,
                const FieldNamesAndFields& ... fieldNamesAndFields
            );

            //- Birth specified elements into the mesh. Fields are expanded
            //  accordingly. New field values can be specified manually by
            //  providing alternating name and field arguments. Values for
            //  other fields will be mapped from the parent elements.
            template<class ... FieldNamesAndFields>
            LagrangianSubMesh birth
            (
                const labelList& parents,
                const FieldNamesAndFields& ... fieldNamesAndFields
            );

            //- Reset the mesh to the old-time conditions
            void reset(const bool initial, const bool final);

            //- Clear all geometry out of the Lagrangian mesh
            void clear();


        // Advanced Modify

            //- Manually partition a sub-set of the mesh given a list of
            //  element-group pairs. Elements not given a group will remain in
            //  a block at the start of the mesh. Elements with multiple groups
            //  will be placed in the block associated with the group with the
            //  largest index. Return offsets to the blocks.
            labelList partition
            (
                const label nGroups,
                const UList<labelPair>& elementsGroups
            );

            //- As above but for enumerated groups
            template<class Enumeration>
            labelList partition
            (
                const label nGroups,
                const UList<elementGroup<Enumeration>>& elementsGroups
            );

            //- Remove specified elements from the mesh. Shuffles everything
            //  else up.
            void remove(const UList<label>& elements);

            //- Remove a specified number of elements from the end of the mesh.
            void remove(const label nElements);


        // Mesh changes

            //- Clear the positions uses during mapping
            virtual void clearPosition();

            //- Store the positions for use during mapping
            virtual void storePosition();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // Write

            //- Write using given format, version and compression
            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write = true
            ) const;

            //- Write settings from the database
            virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LagrangianMeshI.H"

#ifdef NoRepository
    #include "LagrangianMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
