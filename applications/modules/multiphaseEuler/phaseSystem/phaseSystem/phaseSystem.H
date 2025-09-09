/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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
    Class to represent a system of phases.

SourceFiles
    phaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef phaseSystem_H
#define phaseSystem_H

#include "IOdictionary.H"

#include "phaseModel.H"
#include "phaseInterface.H"
#include "phaseInterfaceKey.H"
#include "HashPtrTable.H"
#include "PtrListDictionary.H"
#include "hashedWordList.H"

#include "pimpleNoLoopControl.H"

#include "IOMRFZoneList.H"
#include "fvModels.H"
#include "fvConstraints.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "MULES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class surfaceTensionCoefficientModel;
class momentumTransferSystem;
class pressureReference;
class nonOrthogonalSolutionControl;

/*---------------------------------------------------------------------------*\
                         Class phaseSystem Declaration
\*---------------------------------------------------------------------------*/

class phaseSystem
:
    public IOdictionary
{
public:

    // Public Types

        //- alpha solution control structure
        struct alphaControl
        {
            //- Function to calculate the number of explicit MULES sub-cycles
            //  from the alpha Courant number
            autoPtr<Function1<scalar>> nAlphaSubCyclesPtr;

            //- Number of alpha equation sub-cycles
            label nAlphaSubCycles;

            //- Number of alpha correctors
            //  Usually used to improve the accuracy at high Courant numbers
            //  with semi-implicit MULES, MULESCorr = true
            label nAlphaCorr;

            //- Semi-implicit MULES switch
            bool MULESCorr;

            //- Explicit correction relaxation factor (defaults to 0.5)
            scalar MULESCorrRelax;

            //- Compressibility source stabilisation tolerance
            scalar vDotResidualAlpha;

            //- MULES controls
            MULES::control MULES;

            //- Optionally clip the phase-fractions
            //  between their physical limits and to sum to 1.
            //  Defaults to false
            Switch clip;

            //- Read the alpha and MULES controls from dict
            void read(const dictionary& dict);

            //- Correct nAlphaSubCycles for the current Courant number
            void correct(const scalar CoNum);
        };


    // Public Type Definitions

        //- List of phase models
        typedef PtrListDictionary<phaseModel> phaseModelList;

        //- Partial list of phase models
        typedef UPtrList<phaseModel> phaseModelPartialList;


protected:

    // Protected Type Definitions

        //- Table of interface compression coefficients
        typedef
            HashTable
            <
                scalar,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >
            cAlphaTable;

        //- Table of surface tension models
        typedef
            HashPtrTable
            <
                surfaceTensionCoefficientModel,
                phaseInterfaceKey,
                phaseInterfaceKey::hash
            >
            surfaceTensionCoefficientModelTable;


    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to pimpleNoLoopControl
        const pimpleNoLoopControl& pimple_;

        //- Optional MRF zones
        IOMRFZoneList MRF_;

        //- Name of optional reference phase which is not solved for
        // but obtained from the sum of the other phases
        word referencePhaseName_;

        //- Phase models
        phaseModelList phaseModels_;

        //- Moving phase models
        phaseModelPartialList movingPhaseModels_;

        //- Stationary phase models
        phaseModelPartialList stationaryPhaseModels_;

        //- Thermal phase models
        phaseModelPartialList thermalPhaseModels_;

        //- Multi-component phase models
        phaseModelPartialList multicomponentPhaseModels_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- Rate of change of pressure
        volScalarField dpdt_;

        //- Interface compression coefficients
        const cAlphaTable cAlphas_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- Surface tension models
        const surfaceTensionCoefficientModelTable
            surfaceTensionCoefficientModels_;


    // Protected Member Functions

        //- Return the sum of the phase fractions of the moving phases
        tmp<volScalarField> sumAlphaMoving() const;

        //- Re-normalise the velocity of the phases
        //  around the specified mixture mean
        void setMixtureU(const volVectorField& Um);

        //- Re-normalise the flux of the phases
        //  around the specified mixture mean
        void setMixturePhi
        (
            const PtrList<surfaceScalarField>& alphafs,
            const surfaceScalarField& phim
        );


        // Functions required for interface compression

            //- Normal to interface between two phases
            //  Used for interface compression
            tmp<surfaceVectorField> nHatfv
            (
                const volScalarField& alpha1,
                const volScalarField& alpha2
            ) const;

            //- Normal to interface between two phases dotted with face areas
            //  Used for interface compression
            tmp<surfaceScalarField> nHatf
            (
                const volScalarField& alpha1,
                const volScalarField& alpha2
            ) const;

            //- Curvature of interface between two phases
            //  Used for interface compression
            tmp<volScalarField> K
            (
                const phaseModel& alpha1,
                const phaseModel& alpha2
            ) const;


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

            //- Return pimpleNoLoopControl
            inline const pimpleNoLoopControl& pimple() const;

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
            inline const phaseModelPartialList& thermalPhases() const;

            //- Access the models for phases that have variable temperature
            inline phaseModelPartialList& thermalPhases();

            //- Return the models for phases that have multiple species
            inline const phaseModelPartialList& multicomponentPhases() const;

            //- Access the models for phases that have multiple species
            inline phaseModelPartialList& multicomponentPhases();

            //- Return the phase not given as an argument in a two-phase system
            //  An error is generated if the system is not two-phase
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

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

            //- Access the fvModels
            inline Foam::fvModels& fvModels(fvMesh& mesh);

            //- Access the fvModels
            inline const Foam::fvModels& fvModels() const;

            //- Access the fvConstraints
            inline Foam::fvConstraints& fvConstraints(fvMesh& mesh);

            //- Access the fvConstraints
            inline const Foam::fvConstraints& fvConstraints() const;


        // Sub-model lookup

            //- Check availability of a sub model for a given interface
            template<class ModelType>
            bool foundInterfacialModel
            (
                const phaseInterface& interface
            ) const;

            //- Return a sub model for an interface
            template<class ModelType>
            const ModelType& lookupInterfacialModel
            (
                const phaseInterface& interface
            ) const;


        // Properties

            //- Return the mixture density
            tmp<volScalarField> rho() const;

            //- Return the mixture velocity
            tmp<volVectorField> U() const;

            //- Return the surface tension coefficient for an interface
            tmp<volScalarField> sigma
            (
                const phaseInterfaceKey& key
            ) const;

            //- Return the surface tension coefficient for an interface on a
            //  patch
            tmp<scalarField> sigma
            (
                const phaseInterfaceKey& key,
                const label patchi
            ) const;

            //- Indicator of the proximity of the interface
            //  Field values are 1 near and 0 away for the interface.
            tmp<volScalarField> nearInterface() const;

            //- Stabilisation for normalisation of the interface normal
            inline const dimensionedScalar& deltaN() const;

            //- Return the surface tension force
            tmp<surfaceScalarField> surfaceTension(const phaseModel&) const;

            //- Return incompressibility
            bool incompressible() const;


        // Evolution

            //- Returns true if the phase pressure is treated implicitly
            //  in the phase fraction equations
            bool implicitPhasePressure() const;

            //- Solve for the phase fractions
            void solve
            (
                const alphaControl& alphaControls,
                const PtrList<volScalarField>& rAs,
                const momentumTransferSystem& mts
            );

            //- Correct the fluid properties other than those listed below
            void correct();

            //- Correct the continuity errors
            void correctContinuityError
            (
                const PtrList<volScalarField::Internal>& dmdts
            );

            //- Correct the kinematics
            void correctKinematics();

            //- Correct the thermodynamics
            void correctThermo();

            //- Correct the reactions
            void correctReactions();

            //- Correct the species mass fractions
            void correctSpecies();

            //- Predict the momentumTransport
            void predictMomentumTransport();

            //- Predict the energy transport e.g. alphat
            void predictThermophysicalTransport();

            //- Correct the momentumTransport
            void correctMomentumTransport();

            //- Correct the energy transport e.g. alphat
            void correctThermophysicalTransport();

            //- Update the fluid properties for mesh changes
            void meshUpdate();

            //- Correct fixed-flux BCs to be consistent with the velocity BCs
            void correctBoundaryFlux();

            //- ...
            void correctPhi
            (
                const volScalarField& p_rgh,
                const autoPtr<volScalarField>& divU,
                const pressureReference& pressureReference,
                nonOrthogonalSolutionControl& pimple
            );


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
