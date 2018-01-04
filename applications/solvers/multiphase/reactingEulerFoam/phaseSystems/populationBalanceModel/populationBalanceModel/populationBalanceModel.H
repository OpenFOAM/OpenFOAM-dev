/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::diameterModels::populationBalance

Description
    Class that solves the univariate population balance equation by means of
    a class method (also called sectional or discrete method). The internal
    coordinate is set to the particle volume, so the equation is based on
    a transport equation of the volume-based number density function. The
    discretization is done using the fixed pivot technique of Kumar and
    Ramkrishna (1996). The source terms for coalescence and breakup are written
    in a way that particle number and mass are preserved. An additional breakup
    formulation that is restricted to binary breakup is also included which
    allows a more efficient handling of breakup models with an integrated
    daughter size distribution. It is given by Liao et al. (2017). The
    implementation allows to split the population balance over multiple velocity
    fields using the capability of reactingMultiphaseEulerFoam to solve for n
    momentum equations. It is also possible to define multiple population
    balances, e.g. bubbles and droplets simultaneously.

    References:
    \verbatim
        Coalescence and breakup term formulation:
            Kumar, S. & Ramkrishna, D. (1996).
            On the solution of population balance equations by
            discretization - I. A fixed pivot technique.
            Chemical Engineering Science,
            51(8), 1311-1332.
    \endverbatim

    \verbatim
        Binary breakup term formulation:
            Liao, Y.; Oertel, R.; Kriebitzsch, S.; Schlegel, F.; Lucas, D.
            (2017). A Discrete Population Balance Equation for Binary Breakage.
            International Journal for Numerical Methods in Fluids, In review.
    \endverbatim

Usage
    Example excerpt from a phaseProperties dictionary.
    \verbatim
    type populationBalanceTwoPhaseSystem;

    phases (air water);

    populationBalances (bubbles);

    air
    {
        type            purePhaseModel;
        diameterModel   velocityGroup;
        velocityGroupCoeffs
        {
            populationBalance    bubbles;

            formFactor      0.5235987756;

            sizeGroups
            (
                f0{d  1.00e-3; value 0;}
                f1{d  1.08e-3; value 0;}
                f2{d  1.16e-3; value 0.25;}
                f3{d  1.25e-3; value 0.5;}
                f4{d  1.36e-3; value 0.25;}
                f5{d  1.46e-3; value 0;}
                ...
            );
        }

        residualAlpha   1e-6;
    }

    populationBalanceCoeffs
    {
        bubbles
        {
            continuousPhase water;

            coalescenceModels
            (
                hydrodynamic
                {
                    C 0.25;
                }
            );

            binaryBreakupModels
            ();

            breakupModels
            (
                exponential
                {
                    C 0.5;
                    exponent 0.01;
                    daughterSizeDistributionModel uniformBinary;
                }
            );

            driftModels
            (
                densityChange{}
            );
        }
    }
    \endverbatim

SourceFiles
    populationBalanceModel.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceModel_H
#define populationBalanceModel_H

#include "sizeGroup.H"
#include "phasePair.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseSystem;

namespace diameterModels
{

class coalescenceModel;
class breakupModel;
class binaryBreakupModel;
class driftModel;

/*---------------------------------------------------------------------------*\
                   Class populationBalanceModel Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceModel
:
    public regIOobject
{
    // Private typedefs

        typedef
            HashTable<autoPtr<phasePair>, phasePairKey, phasePairKey::hash>
            phasePairTable;


    // Private data

        //- Reference to the phaseSystem
        const phaseSystem& fluid_;

        //- Interfacial Mass transfer rate between velocityGroups
        HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>& pDmdt_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of the populationBalance
        word name_;

        //- Dictionary
        dictionary dict_;

        //- Reference to pimpleControl
        const pimpleControl& pimple_;

        //- Continuous phase
        const phaseModel& continuousPhase_;

        //- velocityGroups belonging to this populationBalance
        List<velocityGroup*> velocityGroups_;

        //- sizeGroups belonging to this populationBalance
        List<sizeGroup*> sizeGroups_;

        //- List of unordered phasePairs in this populationBalance
        phasePairTable phasePairs_;

        //- sizeGroup boundaries
        PtrList<dimensionedScalar> v_;

        //- Section width required for binary breakup formulation
        PtrList<PtrList<dimensionedScalar>> delta_;

        //- Explicitly treated sources
        PtrList<volScalarField> Su_;

        //- Sources treated implicitly or explicitly depending on sign
        PtrList<volScalarField> SuSp_;

        //- Field for caching pDmdt sources
        autoPtr<volScalarField> Sui_;

        //- Coalescence models
        PtrList<coalescenceModel> coalescence_;

        //- Coalescence rate
        autoPtr<volScalarField> coalescenceRate_;

        //- BreakupModels
        PtrList<breakupModel> breakup_;

        //- Breakup rate
        autoPtr<volScalarField> breakupRate_;

        //- Binary breakup models
        PtrList<binaryBreakupModel> binaryBreakup_;

        //- Binary breakup rate
        autoPtr<volScalarField> binaryBreakupRate_;

        //- Drift models
        PtrList<driftModel> drift_;

        //- Drift rate
        autoPtr<volScalarField> driftRate_;

        //- Ratio between sucessive reprasentative volumes
        autoPtr<volScalarField> r_;

        //- Total void fraction of phases belonging to this populationBalance
        volScalarField alphas_;

        //- Average velocity
        volVectorField U_;

        //- Mean Sauter diameter
        autoPtr<volScalarField> d_;


    // Private member functions

        void registerVelocityGroups();

        void add(sizeGroup* group);

        void createPhasePairs();

        void preSolve();

        dimensionedScalar gamma
        (
            const label i,
            const dimensionedScalar& v
        ) const;

        void birthByCoalescence(const label j, const label k);

        void deathByCoalescence(const label i, const label j);

        void birthByBreakup(const label k, const label model);

        void deathByBreakup(const label i);

        void calcDeltas();

        void birthByBinaryBreakup(const label i, const label j);

        void deathByBinaryBreakup(const label j, const label i);

        void drift(const label i);

        void sources();

        void dmdt();

        void calcAlphas();

        void calcVelocity();

        tmp<volScalarField> dsm() const;


public:

    friend class sizeGroup;
    friend class velocityGroup;

    // Constructor

        populationBalanceModel
        (
            const phaseSystem& fluid,
            const word& name,
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >& pDmdt
        );

        //- Return clone
        autoPtr<populationBalanceModel> clone() const;

        //- Return a pointer to a new populationBalanceModel object created on
        //  freestore from Istream
        class iNew
        {
            const phaseSystem& fluid_;

            HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>&
                pDmdt_;

        public:

            iNew
            (
                const phaseSystem& fluid,
                HashPtrTable<volScalarField, phasePairKey, phasePairKey::hash>&
                    pDmdt
            )
            :
                fluid_(fluid),
                pDmdt_(pDmdt)
            {}

            autoPtr<populationBalanceModel> operator()(Istream& is) const
            {
                return autoPtr<populationBalanceModel>
                (
                    new populationBalanceModel(fluid_, word(is), pDmdt_)
                );
            }
        };


    //- Destructor
    virtual ~populationBalanceModel();

    // Member Functions

        //- Dummy write for regIOobject
        bool writeData(Ostream&) const;

        //- Return reference to the phaseSystem
        inline const phaseSystem& fluid() const;

        //- Return reference to the mesh
        inline const fvMesh& mesh() const;

        //- Return populationBalanceCoeffs dictionary
        inline const dictionary& dict() const;

        //- Return continuous phase
        inline const phaseModel& continuousPhase() const;

        //- Return the velocityGroups belonging to this populationBalance
        inline const List<velocityGroup*>& velocityGroups() const;

        //- Return the sizeGroups belonging to this populationBalance
        inline const List<sizeGroup*>& sizeGroups() const;

        //- Return list of unordered phasePairs in this populationBalance
        inline const phasePairTable& phasePairs() const;

        //- Return the sizeGroup boundaries
        inline const PtrList<dimensionedScalar>& v() const;

        //- Return total void of phases belonging to this populationBalance
        inline const volScalarField& alphas() const;

        //- Return average velocity
        inline const volVectorField& U() const;

        //- Solve the population balance equation
        void solve();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalanceModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
