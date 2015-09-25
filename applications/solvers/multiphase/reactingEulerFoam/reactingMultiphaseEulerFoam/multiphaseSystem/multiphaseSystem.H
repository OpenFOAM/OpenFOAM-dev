/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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
    Foam::multiphaseSystem

Description
    Class which solves the volume fraction equations for two phases.

SourceFiles
    multiphaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseSystem_H
#define multiphaseSystem_H

#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dragModel;
class virtualMassModel;

/*---------------------------------------------------------------------------*\
                         Class multiphaseSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseSystem
:
    public phaseSystem
{
    // Private data

        volScalarField alphas_;

        typedef HashTable<scalar, phasePairKey, phasePairKey::hash>
            cAlphaTable;

        cAlphaTable cAlphas_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        //- Conversion factor for degrees into radians
        static const scalar convertToRad;


    // Private member functions

        void calcAlphas();

        void solveAlphas();

        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        void correctContactAngle
        (
            const phaseModel& alpha1,
            const phaseModel& alpha2,
            surfaceVectorField::GeometricBoundaryField& nHatb
        ) const;

        tmp<volScalarField> K
        (
            const phaseModel& alpha1,
            const phaseModel& alpha2
        ) const;


        //- Return the drag coefficient for phase pair
        virtual tmp<volScalarField> Kd(const phasePairKey& key) const = 0;

        //- Return the face drag coefficient for phase pair
        virtual tmp<surfaceScalarField> Kdf(const phasePairKey& key) const = 0;

        //- Return the virtual mass coefficient for phase pair
        virtual tmp<volScalarField> Vm(const phasePairKey& key) const = 0;

        //- Return the face virtual mass coefficient for phase pair
        virtual tmp<surfaceScalarField> Vmf(const phasePairKey& key) const = 0;

        //- Return the turbulent diffusivity for phase pair
        //  Multiplies the phase-fraction gradient
        virtual tmp<volScalarField> D(const phasePairKey& key) const = 0;

        //- Return the interfacial mass flow rate for phase pair
        virtual tmp<volScalarField> dmdt(const phasePairKey& key) const = 0;


public:

    //- Runtime type information
    TypeName("multiphaseSystem");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            multiphaseSystem,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        multiphaseSystem(const fvMesh&);


    //- Destructor
    virtual ~multiphaseSystem();


    // Selectors

        static autoPtr<multiphaseSystem> New
        (
            const fvMesh& mesh
        );


    // Member Functions

        //- Return the drag coefficient for all phase-pairs
        virtual const phaseSystem::KdTable& Kds() const = 0;

        //- Return the drag coefficient for phase
        virtual tmp<volScalarField> Kd(const phaseModel& phase) const = 0;

        //- Return the combined force (lift + wall-lubrication) for phase pair
        virtual autoPtr<PtrList<Foam::volVectorField> > Fs() const = 0;

        //- Return the turbulent dispersion force on faces for phase pair
        virtual autoPtr<PtrList<Foam::surfaceScalarField> > phiDs
        (
            const PtrList<volScalarField>& rAUs
        ) const = 0;

        //- Return true if there is mass transfer for phase
        virtual bool transfersMass(const phaseModel& phase) const = 0;

        //- Return the total interfacial mass transfer rate for phase
        virtual tmp<volScalarField> dmdt(const phaseModel& phase) const = 0;

        //- Return the momentum transfer matrices
        virtual autoPtr<momentumTransferTable> momentumTransfer() const = 0;

        //- Return the heat transfer matrices
        virtual autoPtr<heatTransferTable> heatTransfer() const = 0;

        //- Return the mass transfer matrices
        virtual autoPtr<massTransferTable> massTransfer() const = 0;

        tmp<surfaceScalarField> surfaceTension(const phaseModel& phase) const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- Solve for the phase fractions
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiphaseSystemI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
