/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Class which models interfacial momenum transfer between a number of phases.
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

template <class modelType>
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
protected:

    // Protected typedefs

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<dragModel> >,
            phasePairKey,
            phasePairKey::hash
        > dragModelTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<virtualMassModel> >,
            phasePairKey,
            phasePairKey::hash
        > virtualMassModelTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<liftModel> >,
            phasePairKey,
            phasePairKey::hash
        > liftModelTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<wallLubricationModel> >,
            phasePairKey,
            phasePairKey::hash
        > wallLubricationModelTable;

        typedef HashTable
        <
            autoPtr<BlendedInterfacialModel<turbulentDispersionModel> >,
            phasePairKey,
            phasePairKey::hash
        > turbulentDispersionModelTable;


private:

    // Private data

        //- Drag coefficients
        phaseSystem::KdTable Kds_;

        //- Virtual mass coefficients
        phaseSystem::VmTable Vms_;

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

            //- Construct element phasei of Fs if not set and return
            //  Used by Fs()
            volVectorField& setF
            (
                PtrList<volVectorField>& Fs, const label phasei
            ) const;

            //- Construct element phasei of phiDs if not set and return
            //  Used by phiDs()
            surfaceScalarField& setPhiD
            (
                PtrList<surfaceScalarField>& phiDs, const label phasei
            ) const;


public:

    // Constructors

        //- Construct from fvMesh
        MomentumTransferPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~MomentumTransferPhaseSystem();


    // Member Functions

        //- Constant access to drag coefficients
        virtual const phaseSystem::KdTable& Kds() const
        {
            return Kds_;
        }

        //- Return the drag coefficient
        virtual tmp<volScalarField> Kd(const phasePairKey& key) const;

        //- Return the face drag coefficient
        virtual tmp<surfaceScalarField> Kdf(const phasePairKey& key) const;

        //- Return the drag coefficient for phase
        virtual tmp<volScalarField> Kd(const phaseModel& phase) const;

        //- Return the virtual mass coefficient
        virtual tmp<volScalarField> Vm(const phasePairKey& key) const;

        //- Return the face virtual mass coefficient
        virtual tmp<surfaceScalarField> Vmf(const phasePairKey& key) const;

        //- Return the combined force (lift + wall-lubrication)
        virtual tmp<volVectorField> F(const phasePairKey& key) const;

        //- Return the combined force (lift + wall-lubrication)
        virtual autoPtr<PtrList<volVectorField> > Fs() const;

        //- Return the turbulent dispersion force on faces for phase pair
        virtual autoPtr<PtrList<surfaceScalarField> > phiDs
        (
            const PtrList<volScalarField>& rAUs
        ) const;

        //- Return the combined face-force (lift + wall-lubrication)
        virtual tmp<surfaceScalarField> Ff(const phasePairKey& key) const;

        //- Return the turbulent diffusivity
        //  Multiplies the phase-fraction gradient
        virtual tmp<volScalarField> D(const phasePairKey& key) const;

        //- Return the momentum transfer matrices
        virtual autoPtr<phaseSystem::momentumTransferTable>
            momentumTransfer() const;

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MomentumTransferPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
