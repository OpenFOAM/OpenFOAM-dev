/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::Reaction

Description
    Simple extension of ReactionThermo to handle reaction kinetics in addition
    to the equilibrium thermodynamics already handled.

SourceFiles
    ReactionI.H
    Reaction.C

\*---------------------------------------------------------------------------*/

#ifndef Reaction_H
#define Reaction_H

#include "reaction.H"
#include "HashPtrTable.H"
#include "scalarField.H"
#include "simpleMatrix.H"
#include "Tuple2.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class ReactionThermo>
class Reaction;

template<class ReactionThermo>
inline Ostream& operator<<(Ostream&, const Reaction<ReactionThermo>&);

class objectRegistry;


/*---------------------------------------------------------------------------*\
                          Class Reaction Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class Reaction
:
    public reaction,
    public ReactionThermo::thermoType
{

public:

    // Static data

        //- Default temperature limits of applicability of reaction rates
        static scalar TlowDefault, ThighDefault;


private:

    // Private Data

        //- Temperature limits of applicability of reaction rates
        scalar Tlow_, Thigh_;


    // Private Member Functions

        //- Construct reaction thermo
        void setThermo(const HashPtrTable<ReactionThermo>& thermoDatabase);


public:

    //- Runtime type information
    TypeName("Reaction");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction,
            dictionary,
            (
                const speciesTable& species,
                const HashPtrTable<ReactionThermo>& thermoDatabase,
                const dictionary& dict
            ),
            (species, thermoDatabase, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction,
            objectRegistry,
            (
                const speciesTable& species,
                const HashPtrTable<ReactionThermo>& thermoDatabase,
                const objectRegistry& ob,
                const dictionary& dict
            ),
            (species, thermoDatabase, ob, dict)
        );


    // Constructors

        //- Construct from components
        Reaction
        (
            const speciesTable& species,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs,
            const HashPtrTable<ReactionThermo>& thermoDatabase
        );

        //- Construct as copy given new speciesTable
        Reaction(const Reaction<ReactionThermo>&, const speciesTable& species);

        //- Construct from dictionary
        Reaction
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ReactionThermo>> clone() const = 0;

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ReactionThermo>> clone
        (
            const speciesTable& species
        ) const = 0;


    // Selectors

        //- Return a pointer to new patchField created on freestore
        //  from dictionary
        static autoPtr<Reaction<ReactionThermo>> New
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const dictionary& dict
        );

        //- Return a pointer to new patchField created on freestore
        //  from objectRegistry and dictionary
        static autoPtr<Reaction<ReactionThermo>> New
        (
            const speciesTable& species,
            const HashPtrTable<ReactionThermo>& thermoDatabase,
            const objectRegistry& ob,
            const dictionary& dict
        );

        //- Return a pointer to new patchField created on freestore from dict.
        //  This is used to construct a single reaction given a speciesThermo
        //  list from which a temporary thermoDatabase is created.
        static autoPtr<Reaction<ReactionThermo>> New
        (
            const speciesTable& species,
            const PtrList<ReactionThermo>& speciesThermo,
            const dictionary& dict
        );


    //- Destructor
    virtual ~Reaction()
    {}


    // Member Functions

        // Access

            //- Return the name of the reaction
            using reaction::name;

            //- Return the lower temperature limit for the reaction
            inline scalar Tlow() const;

            //- Return the upper temperature limit for the reaction
            inline scalar Thigh() const;


        // Hooks

            //- Pre-evaluation hook
            virtual void preEvaluate() const = 0;

            //- Post-evaluation hook
            virtual void postEvaluate() const = 0;


        // Reaction rate coefficients

            //- Net reaction rate for individual species
            void omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const;

            //- Net reaction rate
            scalar omega
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalar& omegaf,
                scalar& omegar
            ) const;

        // Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;

            //- Reverse rate constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;


        // Jacobian coefficients

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the species concentration
            void dwdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarSquareMatrix& J,
                scalarField& dcdt,
                scalar& omegaI,
                scalar& kfwd,
                scalar& kbwd,
                const bool reduced,
                const List<label>& c2s
            ) const;

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the temperature
            void dwdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalar omegaI,
                const scalar kfwd,
                const scalar kbwd,
                scalarSquareMatrix& J,
                const bool reduced,
                const List<label>& c2s,
                const label indexT
            ) const;

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;

            //- Temperature derivative of reverse rate
            virtual scalar dkrdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalar dkfdT,
                const scalar kr
            ) const = 0;

            //- Third-body efficiencies (beta = 1-alpha)
            //  non-empty only for third-body reactions
            //  with enhanced molecularity (alpha != 1)
            virtual const List<Tuple2<label, scalar>>& beta() const = 0;

            //- Species concentration derivative of the pressure dependent term
            virtual void dcidc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dcidc
            ) const = 0;

            //- Temperature derivative of the pressure dependent term
            virtual scalar dcidT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const = 0;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Reaction<ReactionThermo>&) = delete;


    // Ostream Operator

        friend Ostream& operator<< <ReactionThermo>
        (
            Ostream&,
            const Reaction<ReactionThermo>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Reaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
