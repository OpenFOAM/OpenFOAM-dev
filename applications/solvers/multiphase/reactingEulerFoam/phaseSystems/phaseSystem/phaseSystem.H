/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::phaseSystem

Description
    Class to represent a system of phases and model interfacial transfers
    between them.

SourceFiles
    phaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

#include "IOdictionary.H"

#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "HashPtrTable.H"
#include "PtrListDictionary.H"

#include "IOMRFZoneList.H"
#include "fvOptions.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blendingMethod;
template<class modelType> class BlendedInterfacialModel;
class surfaceTensionModel;
class aspectRatioModel;

/*---------------------------------------------------------------------------*\
                         Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary
{
public:

    // Public Typedefs

        typedef HashPtrTable<fvVectorMatrix> momentumTransferTable;

        typedef HashPtrTable<fvScalarMatrix> heatTransferTable;

        typedef HashPtrTable<fvScalarMatrix> specieTransferTable;

        typedef PtrListDictionary<phaseModel> phaseModelList;

        typedef UPtrList<phaseModel> phaseModelPartialList;

        typedef
            HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;

        typedef
            HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>
            dmdtfTable;

        typedef
            HashPtrTable
            <
                HashPtrTable<volScalarField>,
                phasePairKey,
                phasePairKey::hash
            >
            dmidtfTable;



protected:

    // Protected typedefs

        typedef
            HashTable<dictionary, phasePairKey, phasePairKey::hash>
            dictTable;

        typedef
            HashTable<autoPtr<blendingMethod>, word, word::hash>
            blendingMethodTable;

        typedef
            HashTable
            <
                autoPtr<surfaceTensionModel>,
                phasePairKey,
                phasePairKey::hash
            >
            surfaceTensionModelTable;

        typedef
            HashTable
            <
                autoPtr<aspectRatioModel>,
                phasePairKey,
                phasePairKey::hash
            >
            aspectRatioModelTable;


    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Moving phase models
        phaseModelPartialList movingPhaseModels_;

        //- Stationary phase models
        phaseModelPartialList stationaryPhaseModels_;

        //- Anisothermal phase models
        phaseModelPartialList anisothermalPhaseModels_;

        //- Multi-component phase models
        phaseModelPartialList multiComponentPhaseModels_;

        //- Phase pairs
        phasePairTable phasePairs_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- Rate of change of pressure
        volScalarField dpdt_;

        //- Optional MRF zones
        IOMRFZoneList MRF_;

        //- Blending methods
        blendingMethodTable blendingMethods_;


        // Sub Models

            //- Surface tension models
            surfaceTensionModelTable surfaceTensionModels_;

            //- Aspect ratio models
            aspectRatioModelTable aspectRatioModels_;


    // Protected member functions

        //- Calculate and return the mixture flux
        tmp<surfaceScalarField> calcPhi
        (
            const phaseModelList& phaseModels
        ) const;

        //- Generate pairs
        void generatePairs
        (
            const dictTable& modelDicts
        );

        //- Generate pairs and sub-model tables
        template<class modelType>
        void createSubModels
        (
            const dictTable& modelDicts,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models
        );

        //- Generate pairs and sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                autoPtr<modelType>,
                phasePairKey,
                phasePairKey::hash
            >& models,
            const bool correctFixedFluxBCs = true
        );

        //- Generate pairs and blended sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                autoPtr<BlendedInterfacialModel<modelType>>,
                phasePairKey,
                phasePairKey::hash
            >& models,
            const bool correctFixedFluxBCs = true
        );

        //- Generate pairs and two-sided sub-model tables
        template<class modelType>
        void generatePairsAndSubModels
        (
            const word& modelName,
            HashTable
            <
                Pair<autoPtr<modelType>>,
                phasePairKey,
                phasePairKey::hash
            >& models,
            const bool correctFixedFluxBCs = true
        );


public:

    //- Runtime type information
    TypeName("phaseSystem");

    //- Default name of the phase properties dictionary
    static const word propertiesName;


    // Constructors

        //- Construct from fvMesh
        phaseSystem(const fvMesh& mesh);


    //- Destructor
    virtual ~phaseSystem();


    // Member Functions

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the phase models
            inline const phaseModelList& phases() const;

            //- Access the phase models
            inline phaseModelList& phases();

            //- Return the models for phases that are moving
            inline const phaseModelPartialList& movingPhases() const;

            //- Access the models for phases that are moving
            inline phaseModelPartialList& movingPhases();

            //- Return the models for phases that are stationary
            inline const phaseModelPartialList& stationaryPhases() const;

            //- Access the models for phases that are stationary
            inline phaseModelPartialList& stationaryPhases();

            //- Return the models for phases that have variable temperature
            inline const phaseModelPartialList& anisothermalPhases() const;

            //- Access the models for phases that have variable temperature
            inline phaseModelPartialList& anisothermalPhases();

            //- Return the models for phases that have multiple species
            inline const phaseModelPartialList& multiComponentPhases() const;

            //- Access the models for phases that have multiple species
            inline phaseModelPartialList& multiComponentPhases();

            //- Return the phase pairs
            inline const phasePairTable& phasePairs() const;

            //- Return the mixture flux
            inline const surfaceScalarField& phi() const;

            //- Access the mixture flux
            inline surfaceScalarField& phi();

            //- Return the rate of change of the pressure
            inline const volScalarField& dpdt() const;

            //- Access the rate of change of the pressure
            inline volScalarField& dpdt();

            //- Return MRF zones
            inline const IOMRFZoneList& MRF() const;

            //- Access the fvOptions
            inline fv::options& fvOptions() const;


        // Sub-model lookup

            //- Check availability of a sub model for a given phase pair
            template<class modelType>
            bool foundSubModel(const phasePair& key) const;

            //- Return a sub model between a phase pair
            template<class modelType>
            const modelType& lookupSubModel(const phasePair& key) const;

            //- Check availability of a sub model between two phases
            template<class modelType>
            bool foundSubModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous
            ) const;

            //- Return a sub model between two phases
            template<class modelType>
            const modelType& lookupSubModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous
            ) const;

            //- Check availability of a blended sub model for a given phase pair
            template<class modelType>
            bool foundBlendedSubModel(const phasePair& key) const;

            //- Return a blended sub model between a phase pair
            template<class modelType>
            const BlendedInterfacialModel<modelType>&
            lookupBlendedSubModel(const phasePair& key) const;


        // Field construction

            //- Fill up gaps in a phase-indexed list of fields with zeros
            template
            <
                class Type,
                template<class> class PatchField,
                class GeoMesh
            >
            void fillFields
            (
                const word& name,
                const dimensionSet& dims,
                PtrList<GeometricField<Type, PatchField, GeoMesh>>& fieldList
            ) const;

            //- Fill up gaps in a phase-indexed table of fields with zeros
            template
            <
                class Type,
                template<class> class PatchField,
                class GeoMesh
            >
            void fillFields
            (
                const word& name,
                const dimensionSet& dims,
                HashPtrTable<GeometricField<Type, PatchField, GeoMesh>>&
                    fieldTable
            ) const;


        // Properties

            //- Return the mixture density
            tmp<volScalarField> rho() const;

            //- Return the mixture velocity
            tmp<volVectorField> U() const;

            //- Return the aspect-ratio for a pair
            tmp<volScalarField> E(const phasePairKey& key) const;

            //- Return the surface tension coefficient for a pair
            tmp<volScalarField> sigma(const phasePairKey& key) const;

            //- Return the surface tension coefficient for a pair on a patch
            tmp<scalarField> sigma
            (
                const phasePairKey& key,
                label patchi
            ) const;

            //- Return the mass transfer rate for an interface
            virtual tmp<volScalarField> dmdtf(const phasePairKey& key) const;

            //- Return the mass transfer rates for each phase
            virtual PtrList<volScalarField> dmdts() const;


        // Transfers

            //- Return the momentum transfer matrices for the cell-based
            //  algorithm
            virtual autoPtr<momentumTransferTable> momentumTransfer() = 0;

            //- Return the momentum transfer matrices for the face-based
            //  algorithm
            virtual autoPtr<momentumTransferTable> momentumTransferf() = 0;

            //- Return the implicit force coefficients for the face-based
            //  algorithm
            virtual PtrList<surfaceScalarField> AFfs() const = 0;

            //- Return the force fluxes for the cell-based algorithm
            virtual PtrList<surfaceScalarField> phiFs
            (
                const PtrList<volScalarField>& rAUs
            ) = 0;

            //- Return the force fluxes for the face-based algorithm
            virtual PtrList<surfaceScalarField> phiFfs
            (
                const PtrList<surfaceScalarField>& rAUfs
            ) = 0;

            //- Return the force fluxes for the cell-based algorithm
            virtual PtrList<surfaceScalarField> phiKdPhis
            (
                const PtrList<volScalarField>& rAUs
            ) const = 0;

            //- Return the force fluxes for the face-based algorithm
            virtual PtrList<surfaceScalarField> phiKdPhifs
            (
                const PtrList<surfaceScalarField>& rAUfs
            ) const = 0;

            //- Return the explicit part of the drag force
            virtual PtrList<volVectorField> KdUByAs
            (
                const PtrList<volScalarField>& rAUs
            ) const = 0;

            //- Returns true if the phase pressure is treated implicitly
            //  in the phase fraction equation
            virtual bool implicitPhasePressure(const phaseModel& phase) const;

            //- Returns true if the phase pressure is treated implicitly
            //  in the phase fraction equation for any phase
            virtual bool implicitPhasePressure() const;

            //- Return the phase diffusivity
            //  divided by the momentum central coefficient
            virtual PtrList<surfaceScalarField> DByAfs
            (
                const PtrList<volScalarField>& rAUs,
                const PtrList<surfaceScalarField>& rAUfs
            ) const = 0;

            //- Solve the drag system for the new velocities and fluxes
            virtual void partialElimination
            (
                const PtrList<volScalarField>& rAUs
            ) = 0;

            //- Solve the drag system for the new fluxes
            virtual void partialEliminationf
            (
                const PtrList<surfaceScalarField>& rAUfs
            ) = 0;

            //- Return the flux corrections for the cell-based algorithm
            virtual PtrList<surfaceScalarField> ddtCorrByAs
            (
                const PtrList<volScalarField>& rAUs,
                const bool includeVirtualMass = false
            ) const = 0;

            //- Return the heat transfer matrices
            virtual autoPtr<heatTransferTable> heatTransfer() const = 0;

            //- Return the specie transfer matrices
            virtual autoPtr<specieTransferTable> specieTransfer() const = 0;


        // Evolution

            //- Solve for the phase fractions
            virtual void solve
            (
                const PtrList<volScalarField>& rAUs,
                const PtrList<surfaceScalarField>& rAUfs
            );

            //- Correct the fluid properties other than those listed below
            virtual void correct();

            //- Correct the kinematics
            virtual void correctKinematics();

            //- Correct the thermodynamics
            virtual void correctThermo();

            //- Correct the reactions
            virtual void correctReactions();

            //- Correct the species mass fractions
            virtual void correctSpecies();

            //- Correct the turbulence
            virtual void correctTurbulence();

            //- Correct the energy transport e.g. alphat
            virtual void correctEnergyTransport();


        // IO

            //- Read base phaseProperties dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> byDt(const volScalarField& vf);
tmp<surfaceScalarField> byDt(const surfaceScalarField& sf);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
