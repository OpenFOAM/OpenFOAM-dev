/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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
    Foam::solidChemistryModel

Description
    Extends base solid chemistry model by adding a thermo package, and ODE
    functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    solidChemistryModelI.H
    solidChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidChemistryModel_H
#define solidChemistryModel_H

#include "Reaction.H"
#include "ODESystem.H"
#include "volFields.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                   Class solidChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class SolidThermo>
class solidChemistryModel
:
    public CompType,
    public ODESystem
{
    // Private Member Functions

        //- Disallow copy constructor
        solidChemistryModel(const solidChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidChemistryModel&);


protected:

        //- Reference to solid mass fractions
        PtrList<volScalarField>& Ys_;

        //- Reactions
        const PtrList<Reaction<SolidThermo>>& reactions_;

        //- Thermodynamic data of solids
        const PtrList<SolidThermo>& solidThermo_;

        //- Number of solid components
        label nSolids_;

        //- Number of solid reactions
        label nReaction_;

        //- List of reaction rate per solid [kg/m3/s]
        PtrList<volScalarField::Internal> RRs_;

        //- List of active reacting cells
        List<bool> reactingCells_;


    // Protected Member Functions

        //- Write access to source terms for solids
        inline PtrList<volScalarField::Internal>& RRs();

        //- Set reacting status of cell, celli
        void setCellReacting(const label celli, const bool active);


public:

    //- Runtime type information
    TypeName("solidChemistryModel");


    // Constructors

        //- Construct from thermo
        solidChemistryModel(typename CompType::reactionThermo& thermo);


    //- Destructor
    virtual ~solidChemistryModel();


    // Member Functions

        //- The reactions
        inline const PtrList<Reaction<SolidThermo>>& reactions() const;

        //- The number of reactions
        inline label nReaction() const;


        //- dc/dt = omega, rate of change in concentration, for each species
        virtual scalarField omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p,
            const bool updateC0 = false
        ) const = 0;

        //- Return the reaction rate for reaction r and the reference
        //  species and charateristic times
        virtual scalar omega
        (
            const Reaction<SolidThermo>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const = 0;


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
        ) const = 0;

        //- Calculates the reaction rates
        virtual void calculate() = 0;


        // Solid Chemistry model functions

            //- Return const access to the chemical source terms for solids
            inline const volScalarField::Internal& RRs
            (
                const label i
            ) const;

            //- Return total solid source term
            inline tmp<volScalarField::Internal> RRs() const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT) = 0;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s3]
            virtual tmp<volScalarField> Qdot() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            virtual label nEqns() const = 0;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt
            ) const = 0;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                scalarField& dcdt,
                scalarSquareMatrix& dfdc
            ) const = 0;

            virtual void solve
            (
                scalarField &c,
                scalar& T,
                scalar& p,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "solidChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solidChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
