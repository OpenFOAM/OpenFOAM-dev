/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::chemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    chemistryModelI.H
    chemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryModel_H
#define chemistryModel_H

#include "Reaction.H"
#include "ODESystem.H"
#include "volFieldsFwd.H"
#include "simpleMatrix.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class chemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemistryModel
:
    public CompType,
    public ODESystem
{
    // Private Member Functions

        //- Disallow copy constructor
        chemistryModel(const chemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const chemistryModel&);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


protected:

    typedef ThermoType thermoType;

    // Private data

        //- Reference to the field of specie mass fractions
        PtrList<volScalarField>& Y_;

        //- Reactions
        const PtrList<Reaction<ThermoType> >& reactions_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermo_;

        //- Number of species
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- List of reaction rate per specie [kg/m3/s]
        PtrList<DimensionedField<scalar, volMesh> > RR_;


    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<DimensionedField<scalar, volMesh> >& RR();


public:

    //- Runtime type information
    TypeName("chemistryModel");


    // Constructors

        //- Construct from mesh
        chemistryModel(const fvMesh& mesh);


    //- Destructor
    virtual ~chemistryModel();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction<ThermoType> >& reactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<ThermoType>& specieThermo() const;

        //- The number of species
        inline label nSpecie() const;

        //- The number of reactions
        inline label nReaction() const;

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual tmp<scalarField> omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;


        //- Return the reaction rate for iReaction and the reference
        //  species and charateristic times
        virtual scalar omegaI
        (
            label iReaction,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        //- Calculates the reaction rates
        virtual void calculate();


        // Chemistry model functions (overriding abstract functions in
        // basicChemistryModel.H)

            //- Return const access to the chemical source terms for specie, i
            inline const DimensionedField<scalar, volMesh>& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m3/s]
            virtual DimensionedField<scalar, volMesh>& RR
            (
                const label i
            );

            //- Return reaction rate of the specieI in reactionI
            virtual tmp<DimensionedField<scalar, volMesh> > calculateRR
            (
                const label reactionI,
                const label specieI
            ) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return source for enthalpy equation [kg/m/s3]
            virtual tmp<volScalarField> Sh() const;

            //- Return the heat release, i.e. enthalpy/sec [kg/m2/s3]
            virtual tmp<volScalarField> dQ() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const;

            virtual void solve
            (
                scalarField &c,
                scalar& T,
                scalar& p,
                scalar& deltaT,
                scalar& subDeltaT
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
