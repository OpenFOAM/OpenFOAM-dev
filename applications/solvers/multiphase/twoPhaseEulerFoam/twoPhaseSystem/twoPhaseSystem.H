/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::twoPhaseSystem

Description

SourceFiles
    twoPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseSystem_H
#define twoPhaseSystem_H

#include "IOdictionary.H"
#include "phaseModel.H"
#include "phasePair.H"
#include "orderedPhasePair.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class virtualMassModel;
class heatTransferModel;
class liftModel;
class wallLubricationModel;
class turbulentDispersionModel;

class blendingMethod;
template<class modelType> class BlendedInterfacialModel;

/*---------------------------------------------------------------------------*\
                      Class twoPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseSystem
:
    public IOdictionary
{
    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Phase model 1
        phaseModel phase1_;

        //- Phase model 2
        phaseModel phase2_;

        //- Total volumetric flux
        surfaceScalarField phi_;

        //- Dilatation term
        volScalarField dgdt_;

        //- Optional dispersion diffusivity
        tmp<surfaceScalarField> pPrimeByA_;

        //- Unordered phase pair
        autoPtr<phasePair> pair_;

        //- Phase pair for phase 1 dispersed in phase 2
        autoPtr<orderedPhasePair> pair1In2_;

        //- Phase pair for phase 2 dispersed in phase 1
        autoPtr<orderedPhasePair> pair2In1_;

        //- Blending methods
        HashTable<autoPtr<blendingMethod>, word, word::hash> blendingMethods_;

        //- Drag model
        autoPtr<BlendedInterfacialModel<dragModel>> drag_;

        //- Virtual mass model
        autoPtr<BlendedInterfacialModel<virtualMassModel>> virtualMass_;

        //- Heat transfer model
        autoPtr<BlendedInterfacialModel<heatTransferModel>> heatTransfer_;

        //- Lift model
        autoPtr<BlendedInterfacialModel<liftModel>> lift_;

        //- Wall lubrication model
        autoPtr<BlendedInterfacialModel<wallLubricationModel>>
            wallLubrication_;

        //- Wall lubrication model
        autoPtr<BlendedInterfacialModel<turbulentDispersionModel>>
            turbulentDispersion_;


    // Private member functions

        //- Return the mixture flux
        tmp<surfaceScalarField> calcPhi() const;


public:

    // Constructors

        //- Construct from fvMesh
        twoPhaseSystem(const fvMesh&, const dimensionedVector& g);


    //- Destructor
    virtual ~twoPhaseSystem();


    // Member Functions

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture velocity
        tmp<volVectorField> U() const;

        //- Return the drag coefficient
        tmp<volScalarField> Kd() const;

        //- Return the face drag coefficient
        tmp<surfaceScalarField> Kdf() const;

        //- Return the virtual mass coefficient
        tmp<volScalarField> Vm() const;

        //- Return the face virtual mass coefficient
        tmp<surfaceScalarField> Vmf() const;

        //- Return the heat transfer coefficient
        tmp<volScalarField> Kh() const;

        //- Return the combined force (lift + wall-lubrication)
        tmp<volVectorField> F() const;

        //- Return the combined face-force (lift + wall-lubrication)
        tmp<surfaceScalarField> Ff() const;

        //- Return the turbulent diffusivity
        //  Multiplies the phase-fraction gradient
        tmp<volScalarField> D() const;

        //- Solve for the two-phase-fractions
        void solve();

        //- Correct two-phase properties other than turbulence
        void correct();

        //- Correct two-phase turbulence
        void correctTurbulence();

        //- Read base phaseProperties dictionary
        bool read();

        // Access

            //- Access a sub model between a phase pair
            template<class modelType>
            const modelType& lookupSubModel(const phasePair& key) const;

            //- Access a sub model between two phases
            template<class modelType>
            const modelType& lookupSubModel
            (
                const phaseModel& dispersed,
                const phaseModel& continuous
            ) const;

            //- Return the surface tension coefficient
            const dimensionedScalar& sigma() const;

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return phase model 1
            inline const phaseModel& phase1() const;

            //- Return non-const access to phase model 1
            inline phaseModel& phase1();

            //- Return phase model 2
            inline const phaseModel& phase2() const;

            //- Return non-const access to phase model 2
            inline phaseModel& phase2();

            //- Return the phase not given as an argument
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Return the mixture flux
            inline const surfaceScalarField& phi() const;

            //- Return non-const access to the the mixture flux
            inline surfaceScalarField& phi();

            //- Return the dilatation term
            inline const volScalarField& dgdt() const;

            //- Return non-const access to the dilatation parameter
            inline volScalarField& dgdt();

            //- Return non-const access to the dispersion diffusivity
            inline tmp<surfaceScalarField>& pPrimeByA();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "twoPhaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "twoPhaseSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
