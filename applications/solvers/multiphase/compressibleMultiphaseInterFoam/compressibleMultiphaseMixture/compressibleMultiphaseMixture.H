/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    Foam::compressibleMultiphaseMixture

Description

SourceFiles
    compressibleMultiphaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMultiphaseMixture_H
#define compressibleMultiphaseMixture_H

#include "phaseModel.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class compressibleMultiphaseMixture Declaration
\*---------------------------------------------------------------------------*/

class compressibleMultiphaseMixture
:
    public IOdictionary,
    public viscosity
{
public:

    class interfacePair
    :
        public Pair<word>
    {
    public:

        class hash
        :
            public Hash<interfacePair>
        {
        public:

            hash()
            {}

            label operator()(const interfacePair& key) const
            {
                return word::hash()(key.first()) + word::hash()(key.second());
            }
        };


        // Constructors

            interfacePair()
            {}

            interfacePair(const word& alpha1Name, const word& alpha2Name)
            :
                Pair<word>(alpha1Name, alpha2Name)
            {}

            interfacePair(const phaseModel& alpha1, const phaseModel& alpha2)
            :
                Pair<word>(alpha1.name(), alpha2.name())
            {}


        // Friend Operators

            friend bool operator==
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return
                (
                    ((a.first() == b.first()) && (a.second() == b.second()))
                 || ((a.first() == b.second()) && (a.second() == b.first()))
                );
            }

            friend bool operator!=
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (!(a == b));
            }
    };


private:

    // Private Data

        const fvMesh& mesh_;

        //- Pressure
        volScalarField p_;

        //- Mixture temperature
        volScalarField T_;

        //- Dictionary of phases
        PtrDictionary<phaseModel> phases_;

        volScalarField rho_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

        surfaceScalarField rhoPhi_;

        volScalarField alphas_;

        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            sigmaTable;

        sigmaTable sigmas_;
        dimensionSet dimSigma_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;


    // Private Member Functions

        void calcAlphas();

        void solveAlphas(const scalar cAlpha);

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

        tmp<volScalarField> K
        (
            const phaseModel& alpha1,
            const phaseModel& alpha2
        ) const;


public:

    //- Runtime type information
    TypeName("compressibleMultiphaseMixture");


    // Constructors

        //- Construct from components
        compressibleMultiphaseMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~compressibleMultiphaseMixture()
    {}


    // Member Functions

        //- Return the phases
        const PtrDictionary<phaseModel>& phases() const
        {
            return phases_;
        }

        //- Return non-const access to the phases
        PtrDictionary<phaseModel>& phases()
        {
            return phases_;
        }

        //- Return pressure [Pa]
        volScalarField& p()
        {
            return p_;
        }

        //- Return mixture temperature [K]
        volScalarField& T()
        {
            return T_;
        }

        //- Return mixture density [kg/m^3]
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return mixture velocity
        const volVectorField& U() const
        {
            return U_;
        }

        //- Return the volumetric flux
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        const surfaceScalarField& rhoPhi() const
        {
            return rhoPhi_;
        }

        //- Correct the thermodynamics of each phase
        virtual void correctThermo();

        //- Update properties
        virtual void correct();

        //- Update densities for given pressure change
        void correctRho(const volScalarField& dp);

        //- Kinematic viscosity of mixture [m^2/s]
        virtual tmp<volScalarField> nu() const;

        //- Kinematic viscosity of mixture for patch [m^2/s]
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Effective thermal diffusivity of mixture [W/m/K]
        virtual tmp<volScalarField> alphaEff
        (
            const volScalarField& alphat
        ) const;

        //- Return the phase-averaged reciprocal Cv
        tmp<volScalarField> rCv() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- Solve for the mixture phase-fractions
        void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
