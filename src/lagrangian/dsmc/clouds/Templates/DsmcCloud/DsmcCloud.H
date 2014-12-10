/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::DsmcCloud

Description
    Templated base class for dsmc cloud

SourceFiles
    DsmcCloudI.H
    DsmcCloud.C

\*---------------------------------------------------------------------------*/

#ifndef DsmcCloud_H
#define DsmcCloud_H

#include "Cloud.H"
#include "DsmcBaseCloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class BinaryCollisionModel;

template<class CloudType>
class WallInteractionModel;

template<class CloudType>
class InflowBoundaryModel;

/*---------------------------------------------------------------------------*\
                         Class DsmcCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class DsmcCloud
:
    public Cloud<ParcelType>,
    public DsmcBaseCloud
{
    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word cloudName_;

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- A list of unique instances of molecule types in the
        //  simulation.  The position of an entry in the list maps to
        //  the label identifying the typeId, i.e. where typeIdList_ =
        //  (N2 O2 CO2) N2 has typeId label = 0, O2 = 1, CO2 = 2.
        List<word> typeIdList_;

        //- Number of real atoms/molecules represented by a parcel
        scalar nParticle_;

        //- A data structure holding which particles are in which cell
        List<DynamicList<ParcelType*> > cellOccupancy_;

        //- A field holding the value of (sigmaT * cR)max for each
        //  cell (see Bird p220). Initialised with the parcels,
        //  updated as required, and read in on start/restart.
        volScalarField sigmaTcRMax_;

        //- A field holding the remainder from the previous collision selections
        scalarField collisionSelectionRemainder_;

        //- Heat flux at surface field
        volScalarField q_;

        //- Force density at surface field
        volVectorField fD_;

        //- number density field
        volScalarField rhoN_;

        //- Mass density field
        volScalarField rhoM_;

        //- dsmc particle density field
        volScalarField dsmcRhoN_;

        //- linear kinetic energy density field
        volScalarField linearKE_;

        //- Internal energy density field
        volScalarField internalE_;

        // Internal degree of freedom density field
        volScalarField iDof_;

        //- Momentum density field
        volVectorField momentum_;

        //- Parcel constant properties - one for each type
        List<typename ParcelType::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;


        // boundary value fields

            //- boundary temperature
            volScalarField boundaryT_;

            //- boundary velocity
            volVectorField boundaryU_;


        // References to the cloud sub-models

            //- Binary collision model
            autoPtr<BinaryCollisionModel<DsmcCloud<ParcelType> > >
                binaryCollisionModel_;

            //- Wall interaction model
            autoPtr<WallInteractionModel<DsmcCloud<ParcelType> > >
                wallInteractionModel_;

            //- Inflow boundary model
            autoPtr<InflowBoundaryModel<DsmcCloud<ParcelType> > >
                inflowBoundaryModel_;


    // Private Member Functions

        //- Build the constant properties for all of the species
        void buildConstProps();

        //- Record which particles are in which cell
        void buildCellOccupancy();

        //- Initialise the system
        void initialise(const IOdictionary& dsmcInitialiseDict);

        //- Calculate collisions between molecules
        void collisions();

        //- Reset the data accumulation field values to zero
        void resetFields();

        //- Calculate the volume field data
        void calculateFields();

        //- Disallow default bitwise copy construct
        DsmcCloud(const DsmcCloud&);

        //- Disallow default bitwise assignment
        void operator=(const DsmcCloud&);


public:

    // Constructors

        //- Construct given name and mesh, will read Parcels and fields from
        //  file
        DsmcCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            bool readFields = true
        );

        //- Construct given name, mesh and initialisation dictionary.
        DsmcCloud
        (
            const word& cloudName,
            const fvMesh& mesh,
            const IOdictionary& dsmcInitialiseDict
        );


    //- Destructor
    virtual ~DsmcCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            // References to the mesh and databases

                //- Return the cloud type
                inline const word& cloudName() const;

                //- Return refernce to the mesh
                inline const fvMesh& mesh() const;


            // References to the dsmc specific data

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return the idList
                inline const List<word>& typeIdList() const;

                //- Return the number of real particles represented by one
                //  parcel
                inline scalar nParticle() const;

                //- Return the cell occupancy addressing
                inline const List<DynamicList<ParcelType*> >&
                    cellOccupancy() const;

                //- Return the sigmaTcRMax field.  non-const access to allow
                // updating.
                inline volScalarField& sigmaTcRMax();

                //- Return the collision selection remainder field.  non-const
                // access to allow updating.
                inline scalarField& collisionSelectionRemainder();

                //- Return all of the constant properties
                inline const List<typename ParcelType::constantProperties>&
                    constProps() const;

                //- Return the constant properties of the given typeId
                inline const typename ParcelType::constantProperties&
                    constProps(label typeId) const;

                //- Return refernce to the random object
                inline Random& rndGen();


            // References to the boundary fields for surface data collection

                //- Return non-const heat flux boundary field reference
                inline volScalarField::GeometricBoundaryField& qBF();

                //- Return non-const force density at boundary field reference
                inline volVectorField::GeometricBoundaryField& fDBF();

                //- Return non-const number density boundary field reference
                inline volScalarField::GeometricBoundaryField& rhoNBF();

                //- Return non-const mass density boundary field reference
                inline volScalarField::GeometricBoundaryField& rhoMBF();

                //- Return non-const linear kinetic energy density boundary
                //  field reference
                inline volScalarField::GeometricBoundaryField& linearKEBF();

                //- Return non-const internal energy density boundary field
                // reference
                inline volScalarField::GeometricBoundaryField& internalEBF();

                //- Return non-const internal degree of freedom density boundary
                //  field reference
                inline volScalarField::GeometricBoundaryField& iDofBF();

                //- Return non-const momentum density boundary field reference
                inline volVectorField::GeometricBoundaryField& momentumBF();


            // References to the macroscopic fields

                //- Return macroscopic temperature
                inline const volScalarField& boundaryT() const;

                //- Return macroscopic velocity
                inline const volVectorField& boundaryU() const;

                //- Return heat flux at surface field
                inline const volScalarField& q() const;

                //- Return force density at surface field
                inline const volVectorField& fD() const;

                //- Return the real particle number density field
                inline const volScalarField& rhoN() const;

                //- Return the particle mass density field
                inline const volScalarField& rhoM() const;

                //- Return the field of number of DSMC particles
                inline const volScalarField& dsmcRhoN() const;

                //- Return the total linear kinetic energy (translational and
                // thermal density field
                inline const volScalarField& linearKE() const;

                //- Return the internal energy density field
                inline const volScalarField& internalE() const;

                //- Return the average internal degrees of freedom  field
                inline const volScalarField& iDof() const;

                //- Return the momentum density field
                inline const volVectorField& momentum() const;


            // Kinetic theory helper functions

                //- Generate a random velocity sampled from the Maxwellian speed
                // distribution
                vector equipartitionLinearVelocity
                (
                    scalar temperature,
                    scalar mass
                );

                //- Generate a random internal energy, sampled from the
                // equilibrium distribution (Bird eqn 11.22 and 11.23 and
                // adapting code from DSMC3.FOR)
                scalar equipartitionInternalEnergy
                (
                    scalar temperature,
                    scalar internalDegreesOfFreedom
                );


                // From the Maxwellian distribution:
                //- Average particle speed
                inline scalar maxwellianAverageSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianAverageSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                //- RMS particle speed
                inline scalar maxwellianRMSSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianRMSSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                //- Most probable speed
                inline scalar maxwellianMostProbableSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianMostProbableSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;


            // Sub-models

                //- Return reference to binary elastic collision model
                inline const BinaryCollisionModel<DsmcCloud<ParcelType> >&
                    binaryCollision() const;

                //- Return non-const reference to binary elastic collision model
                inline BinaryCollisionModel<DsmcCloud<ParcelType> >&
                    binaryCollision();

                //- Return reference to wall interaction model
                inline const WallInteractionModel<DsmcCloud<ParcelType> >&
                    wallInteraction() const;

                //- Return non-const reference to wall interaction model
                inline WallInteractionModel<DsmcCloud<ParcelType> >&
                    wallInteraction();

                //- Return reference to wall interaction model
                inline const InflowBoundaryModel<DsmcCloud<ParcelType> >&
                    inflowBoundary() const;

                //- Return non-const reference to wall interaction model
                inline InflowBoundaryModel<DsmcCloud<ParcelType> >&
                    inflowBoundary();


        // Check

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total internal energy in the system
            inline scalar internalEnergyOfSystem() const;

            //- Print cloud information
            void info() const;

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;




        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const vector& U,
                const scalar Ei,
                const label cellI,
                const label tetFaceI,
                const label tetPtI,
                const label typeId
            );

            //- Evolve the cloud (move, collide)
            void evolve();

            //- Clear the Cloud
            inline void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DsmcCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DsmcCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
