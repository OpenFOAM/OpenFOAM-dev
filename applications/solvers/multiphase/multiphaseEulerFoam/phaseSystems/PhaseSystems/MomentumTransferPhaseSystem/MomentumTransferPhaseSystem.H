/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::MomentumTransferPhaseSystem

Description
    Class which models interfacial momentum transfer between a number of phases.
    Drag, virtual mass, lift, wall lubrication and turbulent dispersion are all
    modelled. The explicit contribution from the drag is omitted from the
    transfer matrices, as this forms part of the solution of the pressure
    equation.

SourceFiles
    MomentumTransferPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef MomentumTransferPhaseSystem_H
#define MomentumTransferPhaseSystem_H

#include "phaseSystem.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class modelType>
class BlendedInterfacialModel;

class blendingMethod;
class dragModel;
class virtualMassModel;
class liftModel;
class wallLubricationModel;
class turbulentDispersionModel;

/*---------------------------------------------------------------------------*\
                 Class MomentumTransferPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class MomentumTransferPhaseSystem
:
    public BasePhaseSystem
{
private:

    // Private typedefs

        typedef HashPtrTable
        <
            volScalarField,
            phasePairKey,
            phasePairKey::hash
        > KdTable;

        typedef HashPtrTable
        <
            surfaceScalarField,
            phasePairKey,
            phasePairKey::hash
        > KdfTable;

        typedef HashPtrTable
        <
            volScalarField,
            phasePairKey,
            phasePairKey::hash
        > VmTable;

        typedef HashPtrTable
        <
            surfaceScalarField,
            phasePairKey,
            phasePairKey::hash
        > VmfTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<dragModel>>,
            phasePairKey,
            phasePairKey::hash
        > dragModelTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<virtualMassModel>>,
            phasePairKey,
            phasePairKey::hash
        > virtualMassModelTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<liftModel>>,
            phasePairKey,
            phasePairKey::hash
        > liftModelTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<wallLubricationModel>>,
            phasePairKey,
            phasePairKey::hash
        > wallLubricationModelTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<turbulentDispersionModel>>,
            phasePairKey,
            phasePairKey::hash
        > turbulentDispersionModelTable;


    // Private Data

        //- Drag coefficients
        KdTable Kds_;

        //- Face drag coefficients
        KdfTable Kdfs_;

        //- Virtual mass coefficients
        VmTable Vms_;

        //- Face virtual mass coefficients
        VmfTable Vmfs_;


        // Sub Models

            //- Drag models
            dragModelTable dragModels_;

            //- Virtual mass models
            virtualMassModelTable virtualMassModels_;

            //- Lift models
            liftModelTable liftModels_;

            //- Wall lubrication models
            wallLubricationModelTable wallLubricationModels_;

            //- Turbulent dispersion models
            turbulentDispersionModelTable turbulentDispersionModels_;


protected:

    // Protected Member Functions

        //- Add momentum transfer terms which result from bulk mass transfers
        void addDmdtUfs
        (
            const phaseSystem::dmdtfTable& dmdtfs,
            phaseSystem::momentumTransferTable& eqns
        );

        //- Return the drag coefficient for the phase pair
        virtual tmp<volScalarField> Kd(const phasePairKey& key) const;

        //- Return the face drag coefficient for the phase pair
        virtual tmp<surfaceScalarField> Kdf(const phasePairKey& key) const;

        //- Return the virtual mass coefficient for the phase pair
        virtual tmp<volScalarField> Vm(const phasePairKey& key) const;


public:

    // Constructors

        //- Construct from fvMesh
        MomentumTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~MomentumTransferPhaseSystem();


    // Member Functions

        //- Return the momentum transfer matrices for the cell-based algorithm.
        //  This includes implicit and explicit forces that add into the cell
        //  UEqn in the normal way.
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransfer();

        //- As momentumTransfer, but for the face-based algorithm
        virtual autoPtr<phaseSystem::momentumTransferTable> momentumTransferf();

        //- Return implicit force coefficients on the faces, for the face-based
        //  algorithm.
        virtual PtrList<surfaceScalarField> AFfs() const;

        //- Return the explicit force fluxes for the cell-based algorithm, that
        //  do not depend on phase mass/volume fluxes, and can therefore be
        //  evaluated outside the corrector loop. This includes things like
        //  lift, turbulent dispersion, and wall lubrication.
        virtual PtrList<surfaceScalarField> phiFs
        (
            const PtrList<volScalarField>& rAUs
        );

        //- As phiFs, but for the face-based algorithm
        virtual PtrList<surfaceScalarField> phiFfs
        (
            const PtrList<surfaceScalarField>& rAUfs
        );

        //- Return the explicit drag force fluxes for the cell-based algorithm.
        //  These depend on phase mass/volume fluxes, and must therefore be
        //  evaluated inside the corrector loop.
        virtual PtrList<surfaceScalarField> phiKdPhis
        (
            const PtrList<volScalarField>& rAUs
        ) const;

        //- As phiKdPhis, but for the face-based algorithm
        virtual PtrList<surfaceScalarField> phiKdPhifs
        (
            const PtrList<surfaceScalarField>& rAUfs
        ) const;

        //- Return the explicit part of the drag force for the cell-based
        //  algorithm. This is the cell-equivalent of phiKdPhis. These depend on
        //  phase velocities, and must therefore be evaluated inside the
        //  corrector loop.
        virtual PtrList<volVectorField> KdUByAs
        (
            const PtrList<volScalarField>& rAUs
        ) const;

        //- Solve the drag system for the velocities and fluxes
        virtual void partialElimination
        (
            const PtrList<volScalarField>& rAUs,
            const PtrList<volVectorField>& KdUByAs,
            const PtrList<surfaceScalarField>& alphafs,
            const PtrList<surfaceScalarField>& phiKdPhis
        );

        //- As partialElimination, but for the face-based algorithm. Only solves
        //  for the fluxes.
        virtual void partialEliminationf
        (
            const PtrList<surfaceScalarField>& rAUfs,
            const PtrList<surfaceScalarField>& alphafs,
            const PtrList<surfaceScalarField>& phiKdPhifs
        );

        //- Return the flux corrections for the cell-based algorithm. These
        //  depend on phase mass/volume fluxes, and must therefore be evaluated
        //  inside the corrector loop.
        virtual PtrList<surfaceScalarField> ddtCorrByAs
        (
            const PtrList<volScalarField>& rAUs,
            const bool includeVirtualMass = false
        ) const;

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
        ) const;

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MomentumTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
