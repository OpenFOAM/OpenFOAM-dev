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
    Foam::cloud

Description
    Base class for clouds. Provides a basic evolution algorithm, models, and a
    database for caching derived and average fields to avoid unnecessary
    re-calculation.

SourceFiles
    cloud.C
    cloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cloud_H
#define cloud_H

#include "LagrangianModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class CloudStateField;

template<class Type>
class CloudDerivedField;

template<class Type>
class CloudAverageField;

/*---------------------------------------------------------------------------*\
                           Class cloud Declaration
\*---------------------------------------------------------------------------*/

class cloud
:
    public regIOobject
{
private:

    // Private Data

        //- Reference to the mesh
        LagrangianMesh& mesh_;

        //- The demand-driven Lagrangian models pointer
        mutable Foam::LagrangianModels* LagrangianModelsPtr_;

        //- State fields
        #define DECLARE_STATE_FIELDS(Type, nullArg)       \
            mutable PtrList<CloudStateField<Type>>        \
                CAT3(state, CAPITALIZE(Type), Fields_);
        FOR_ALL_FIELD_TYPES(DECLARE_STATE_FIELDS);
        #undef DECLARE_STATE_FIELDS

        //- Derived fields
        #define DECLARE_DERIVED_FIELDS(Type, nullArg)     \
            mutable PtrList<CloudDerivedField<Type>>      \
                CAT3(derived, CAPITALIZE(Type), Fields_);
        FOR_ALL_FIELD_TYPES(DECLARE_DERIVED_FIELDS);
        #undef DECLARE_DERIVED_FIELDS

        //- Average fields
        #define DECLARE_AVERAGE_FIELDS(Type, nullArg)     \
            mutable PtrList<CloudAverageField<Type>>      \
                CAT3(average, CAPITALIZE(Type), Fields_);
        FOR_ALL_FIELD_TYPES(DECLARE_AVERAGE_FIELDS);
        #undef DECLARE_AVERAGE_FIELDS

        //- Retained state labels
        autoPtr<LagrangianLabelInternalField> statePtr_;

        //- Cell length scale
        scalarField cellLengthScaleVf_;


    // Private Member Functions

        //- Construct and store the Lagrangian mesh, and return a reference
        static LagrangianMesh& mesh
        (
            const polyMesh& pMesh,
            const word& name,
            const IOobject::readOption readOption=IOobject::READ_IF_PRESENT,
            const IOobject::writeOption writeOption=IOobject::AUTO_WRITE
        );

        //- Access the state fields. Only specialisations are defined.
        template<class Type>
        PtrList<CloudStateField<Type>>& stateFields() const;

        //- Clear the state fields
        void clearStateFields();

        //- Access the derived fields. Only specialisations are defined.
        template<class Type>
        PtrList<CloudDerivedField<Type>>& derivedFields() const;

        //- Clear the derived fields
        void clearDerivedFields(const bool final);

        //- Access the average fields. Only specialisations are defined.
        template<class Type>
        PtrList<CloudAverageField<Type>>& averageFields() const;

        //- Remove values from a sub mesh from the average fields
        void removeFromAverageFields(const LagrangianSubMesh&);

        //- Add values from a sub mesh into the average fields
        void addToAverageFields(const LagrangianSubMesh&, const bool final);

        //- Correct values from a sub mesh in the average fields
        void correctAverageFields(const LagrangianSubMesh&, const bool final);

        //- Clear the average fields
        void clearAverageFields();

        //- Reset the average fields
        void resetAverageFields();

        //- Return the IO for the retained state labels
        IOobject stateIo(const IOobject::readOption) const;

        //- Read and return the retained state labels, if the file exists,
        //  otherwise return an empty pointer
        autoPtr<LagrangianLabelInternalField> readStates() const;

        //- Return the retained state enumerations, if they are available,
        //  otherwise return an empty pointer
        autoPtr<List<LagrangianState>> initialStates() const;

        //- Store or clear the retained state labels as necessary
        bool storeStates();

        //- Store the retained state labels for the given sub-mesh
        void storeStates(const LagrangianSubMesh& subMesh);

        //- Calculate the cell length scale for the given sub-mesh
        tmp<LagrangianSubScalarField> cellLengthScale
        (
            const LagrangianSubMesh& subMesh
        ) const;

        //- Track through the mesh
        void track
        (
            LagrangianSubScalarSubField& fraction,
            const scalar maxTimeStepFraction,
            const scalar maxCellLengthScaleFraction
        );

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream&) const;


protected:

    // Protected Member Functions

        //- Initialisation hook
        virtual void initialise(const bool predict);

        //- Partition hook
        virtual void partition();

        //- Return the acceleration with which to do second-order tracking
        virtual tmp<LagrangianSubVectorField> dUdt
        (
            const LagrangianSubMesh&
        ) const = 0;

        //- Do we need to re-calculate particles that are modified?
        virtual bool reCalculateModified() = 0;

        //- Update the cloud properties
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );


public:

    // Public Static Data

        //- Run-time type information
        TypeName("cloud");


    // Public Data

        //- Context in which this cloud is used
        const enum contextType { unknown, functionObject, fvModel } context;

        //- Tracking method
        const enum trackingType { linear, parabolic } tracking;

        //- Velocity
        CloudStateField<vector>& U;


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        cloud,
        polyMesh,
        (const polyMesh& mesh, const word& name, const contextType context),
        (mesh, name, context)
    );


    // Constructors

        //- Construct from a mesh and a name
        cloud
        (
            const polyMesh& mesh,
            const word& name,
            const contextType context,
            const IOobject::readOption readOption=IOobject::READ_IF_PRESENT,
            const IOobject::writeOption writeOption=IOobject::AUTO_WRITE
        );

        //- Disallow default bitwise copy construction
        cloud(const cloud&) = delete;


    //- Selector
    static autoPtr<cloud> New
    (
        const polyMesh& mesh,
        const word& name,
        const word& type
    );


    //- Destructor
    virtual ~cloud();


    // Member Functions

        // Access

            //- Access the mesh
            inline const LagrangianMesh& mesh() const
            {
                return mesh_;
            }

            //- Lookup the cloud associated with a mesh
            static const Foam::cloud& lookup(const LagrangianMesh& mesh);

            //- Access the models
            Foam::LagrangianModels& LagrangianModels() const;


        // Fields

            //- Add/get a state field to/from the cache
            template<class Type, class ... Args>
            CloudStateField<Type>& stateField
            (
                const Args& ... args
            ) const;

            //- Add/get a derived field to/from the cache
            template<class Type, class ... Args>
            const CloudDerivedField<Type>& derivedField
            (
                const Args& ... args
            ) const;

            //- Add/get an average field to/from the cache
            template<class Type, class ... Args>
            const CloudAverageField<Type>& averageField
            (
                const Args& ... args
            ) const;


        // Evolution

            //- Solve the cloud's evolution over the current time-step
            virtual void solve();


        // Mesh changes

            //- Store the positions for use during mapping
            virtual void storePosition();

            //- Update for mesh motion
            virtual void movePoints(const polyMesh&);

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianMesh&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cloudTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
