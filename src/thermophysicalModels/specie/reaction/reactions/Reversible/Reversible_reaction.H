/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::reaction::Reversible

Description
    Extension of reaction to handle reversible reactions expressed as a rate
    constant multiplied by powers of the reactants' concentrations. The forward
    rate constants is specified, and the reverse rate constant is determined
    from equilibrium thermodynamics.

SourceFiles
    Reversible_reaction.C

\*---------------------------------------------------------------------------*/

#ifndef Reversible_reaction_H
#define Reversible_reaction_H

#include "reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactions
{

/*---------------------------------------------------------------------------*\
                     Class Reversible Declaration
\*---------------------------------------------------------------------------*/

template<class Rate, class ThermoType>
class Reversible
:
    public reaction
{
public:

    // Public Data Types

        //- Inherit the state structure
        using typename reaction::state;


private:

    // Private Data

        //- Forward rate constant
        Rate k_;

        //- The reaction thermo
        typename ThermoType::thermoType thermo_;


    // Private Member Functions

        //- Construct reaction thermo
        static typename ThermoType::thermoType thermo
        (
            const reactionCoeffs& r,
            const SpecieThermosRef<ThermoType>& specieThermos
        );


public:

    //- Runtime type information
    TypeName("Reversible");


    // Constructors

        //- Construct from components
        Reversible
        (
            const reactionCoeffs& r,
            const specieThermosRef& specieThermos,
            const Rate& k
        );

        //- Construct from dictionary
        Reversible
        (
            const speciesTable& species,
            const specieThermosRef& specieThermos,
            const dictionary& dict
        );

        //- Construct from objectRegistry and dictionary
        Reversible
        (
            const speciesTable& species,
            const specieThermosRef& specieThermos,
            const objectRegistry& ob,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Reversible()
    {}


    // Member Functions

        // Access

            //- Return the selection type
            virtual word selectionType() const;

            //- Return this reaction's Jacobian density
            virtual reaction::jacobianDensity density() const;


        // Hooks

            //- Pre-evaluation hook
            virtual void preEvaluate() const;

            //- Post-evaluation hook
            virtual void postEvaluate() const;


        // Reaction rates

            //- Net reaction rate
            virtual scalar omega
            (
                const state& s,
                scalar& omegaf,
                scalar& omegar
            ) const;

            //- Net reaction rate and derivatives
            virtual scalar omega
            (
                const state& s,
                scalar& dOmegaDT,
                scalarField& dOmegaDc
            ) const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace reactions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Reversible_reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
