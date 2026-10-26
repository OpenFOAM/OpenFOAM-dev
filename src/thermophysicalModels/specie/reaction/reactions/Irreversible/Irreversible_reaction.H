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
    Foam::reaction::Irreversible

Description
    Extension of reaction to handle irreversible reactions expressed as a rate
    constant multiplied by powers of the reactants' concentrations. The forward
    rate constant is specified.

SourceFiles
    Irreversible_reaction.C

\*---------------------------------------------------------------------------*/

#ifndef Irreversible_reaction_H
#define Irreversible_reaction_H

#include "reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace reactions
{

/*---------------------------------------------------------------------------*\
                    Class Irreversible Declaration
\*---------------------------------------------------------------------------*/

template<class Rate>
class Irreversible
:
    public reaction
{
public:

    // Public Data Types

        //- Inherit the state structure
        using typename reaction::state;


    // Private Data

        //- Forward rate constant
        Rate k_;


public:

    //- Runtime type information
    TypeName("Irreversible");


    // Constructors

        //- Construct from components
        Irreversible
        (
            const reactionCoeffs& r,
            const specieThermosRef& specieThermos,
            const Rate& reactionRate
        );

        //- Construct from dictionary
        Irreversible
        (
            const speciesTable& species,
            const specieThermosRef& specieThermos,
            const dictionary& dict
        );

        //- Construct from objectRegistry and dictionary
        Irreversible
        (
            const speciesTable& species,
            const specieThermosRef& specieThermos,
            const objectRegistry& ob,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Irreversible()
    {}


    // Member Functions

        // Access

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
    #include "Irreversible_reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
