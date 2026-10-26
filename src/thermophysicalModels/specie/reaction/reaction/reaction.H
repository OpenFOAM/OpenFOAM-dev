/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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
    Foam::reaction

Description
    Base-class for reactions

SourceFiles
    reactionI.H
    reaction.C
    reactionNew.C

\*---------------------------------------------------------------------------*/

#ifndef reaction_H
#define reaction_H

#include "reactionCoeffs.H"
#include "specieThermosRef.H"
#include "scalarMatrices.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class reaction Declaration
\*---------------------------------------------------------------------------*/

class reaction
:
    public reactionCoeffs
{
public:

    // Public Data Types

        //- Reaction state structure
        struct state;

        //- A "no" complete-to-simplified map. A hard-coded identity map.
        struct noc2s;

        //- Enumeration of the reactant's Jacobian density
        enum class jacobianDensity
        {
            reactants,
            reactantsAndProducts,
            all
        };


private:

    // Private Member Functions

        //- Map from old/backwards-compatible types to current types
        static HashTable<word> backwardTypeToType();

        //- Convert from an old/backwards-compatible type to a current type
        static word backwardTypeToType(const word&);

        //- Dynamically compile a thermo reaction
        static bool compile
        (
            const specieThermosRef& specieThermos,
            ITstream& typeIts,
            const word& type,
            const word& typeThermo,
            const word& table
        );


protected:

    // Protected Member Data

        //- Work array the size of the left hand side
        mutable scalarField lWork_;

        //- Work array the size of the right hand side
        mutable scalarField rWork_;


    // Protected Member Functions

        //- Dimensions of the rate constant
        static inline dimensionSet kDimensions
        (
            const List<specieCoeffs>& coeffs
        );

        //- Return the product of the concentrations raised to their exponents
        static inline scalar cPowExponentProduct
        (
            const List<specieCoeffs>& coeffs,
            const state& s
        );

        //- Return the product of the concentrations raised to their exponents
        static inline scalar cPowExponentProduct
        (
            const List<specieCoeffs>& coeffs,
            const state& s,
            scalarField& dResultDc
        );


public:

    //- Runtime type information
    TypeName("reaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            reaction,
            dictionary,
            (
                const speciesTable& species,
                const specieThermosRef& specieThermos,
                const dictionary& dict
            ),
            (species, specieThermos, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            reaction,
            objectRegistry,
            (
                const speciesTable& species,
                const specieThermosRef& specieThermos,
                const objectRegistry& ob,
                const dictionary& dict
            ),
            (species, specieThermos, ob, dict)
        );


    // Constructors

        //- Construct from coefficients
        reaction(const reactionCoeffs&);

        //- Construct from dictionary
        reaction
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Selectors

        //- Return a pointer to new reaction created from a dictionary
        static autoPtr<reaction> New
        (
            const speciesTable& species,
            const specieThermosRef& specieThermos,
            const dictionary& dict
        );

        //- Return a pointer to new reaction created from an objectRegistry and
        //  a dictionary
        static autoPtr<reaction> New
        (
            const speciesTable& species,
            const specieThermosRef& specieThermos,
            const objectRegistry& ob,
            const dictionary& dict
        );


    //- Destructor
    virtual ~reaction()
    {}


    // Member Functions

        // Access

            //- Return the selection type
            virtual word selectionType() const;

            //- Return this reaction's Jacobian density
            virtual jacobianDensity density() const = 0;


        // Hooks

            //- Pre-evaluation hook
            virtual void preEvaluate() const = 0;

            //- Post-evaluation hook
            virtual void postEvaluate() const = 0;


        // Reaction rates

            //- Net reaction rate
            inline scalar omega(const state& s) const;

            //- Net reaction rate
            virtual scalar omega
            (
                const state& s,
                scalar& omegaf,
                scalar& omegar
            ) const = 0;

            //- Net reaction rate and derivatives
            virtual scalar omega
            (
                const state& s,
                scalar& dOmegaDT,
                scalarField& dOmegaDc
            ) const = 0;

            //- The net reaction rate for each species involved, given the
            //  total net reaction rate
            template<class C2S>
            inline void dNdtByV
            (
                const scalar omega,
                scalarField& dNdtByV,
                const C2S& c2s
            ) const;

            //- The net reaction rate for each species involved, given the
            //  total net reaction rate
            inline void dNdtByV
            (
                const scalar omega,
                scalarField& dNdtByV,
                const bool reduced,
                const List<label>& c2s
            ) const;

            //- The net reaction rate for each species involved
            template<class C2S>
            inline void dNdtByV
            (
                const state& s,
                scalarField& dNdtByV,
                const C2S& c2s
            ) const;

            //- The net reaction rate for each species involved
            inline void dNdtByV
            (
                const state& s,
                scalarField& dNdtByV,
                const bool reduced,
                const List<label>& c2s
            ) const;


        // Jacobian of reaction rates

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the concentration and temperature
            template<class C2S>
            inline void ddNdtByVdcT
            (
                const state& s,
                scalarField& dNdtByV,
                scalarSquareMatrix& ddNdtByVdcT,
                const C2S& c2s,
                const label Tj
            ) const;

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the concentration and temperature
            virtual inline void ddNdtByVdcT
            (
                const state& s,
                scalarField& dNdtByV,
                scalarSquareMatrix& ddNdtByVdcT,
                const bool reduced,
                const List<label>& c2s,
                const label Tj
            ) const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reaction&) = delete;
};


//- Write to a stream
Ostream& operator<<(Ostream&, const reaction&);


/*---------------------------------------------------------------------------*\
                       Class reaction::state Declaration
\*---------------------------------------------------------------------------*/

struct reaction::state
{
    // Public Member Data

        //- Pressure
        scalar p;

        //- Temperature
        scalar T;

        //- Species' concentrations
        scalarField c;

        //- Index
        label li;

        //- Logarithm of the pressure
        scalar lnp;

        //- Reciprocal of the temperature
        scalar rT;

        //- Logarithm of the temperature
        scalar lnT;

        //- Total concentration of all species
        scalar C;

        //- Work arrays the size of the number of species
        mutable FixedList<scalarField, 3> specieWork_;


    // Constructors

        //- Construct for a number of species
        state(const label nSpecie);


    // Member Functions

        //- Reset all values apart from species' concentrations
        inline void reset
        (
            const scalar p,
            const scalar T,
            const label li
        );

        //- Reset all values
        inline void reset
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li
        );
};


/*---------------------------------------------------------------------------*\
                       Class reaction::state Declaration
\*---------------------------------------------------------------------------*/

struct reaction::noc2s
{
    inline label operator[](const label c) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "reactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
