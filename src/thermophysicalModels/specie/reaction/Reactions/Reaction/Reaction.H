/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Simple extension of ThermoType to handle reaction kinetics in
    addition to the equilibrium thermodynamics already handled.

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

template<class ThermoType>
class Reaction;

template<class ThermoType>
inline Ostream& operator<<(Ostream&, const Reaction<ThermoType>&);

class objectRegistry;


/*---------------------------------------------------------------------------*\
                          Class Reaction Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class Reaction
:
    public reaction,
    public ThermoType::thermoType
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
        void setThermo(const PtrList<ThermoType>& speciesThermo);


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
                const PtrList<ThermoType>& speciesThermo,
                const dictionary& dict
            ),
            (species, speciesThermo, dict)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            Reaction,
            objectRegistry,
            (
                const speciesTable& species,
                const PtrList<ThermoType>& speciesThermo,
                const objectRegistry& ob,
                const dictionary& dict
            ),
            (species, speciesThermo, ob, dict)
        );


    // Constructors

        //- Construct from components
        Reaction
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs
        );

        //- Construct as copy given new speciesTable
        Reaction
        (
            const Reaction<ThermoType>&,
            const speciesTable& species
        );

        //- Construct from dictionary
        Reaction
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ThermoType>> clone() const = 0;

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ThermoType>> clone
        (
            const speciesTable& species
        ) const = 0;


    // Selectors

        //- Return a pointer to new reaction created from a dictionary
        static autoPtr<Reaction<ThermoType>> New
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const dictionary& dict
        );

        //- Return a pointer to new reaction created from an objectRegistry and
        //  a dictionary
        static autoPtr<Reaction<ThermoType>> New
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const objectRegistry& ob,
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

            //- Concentration powers
            void C
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalar& Cf,
                scalar& Cr
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

            //- The net reaction rate for each species involved
            void dNdtByV
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dNdtByV,
                const bool reduced,
                const List<label>& c2s,
                const label Nsi0
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
                const scalar kf,
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

            //- Does this reaction have concentration-dependent rate constants?
            virtual bool hasDkdc() const = 0;

            //- Concentration derivative of forward rate
            virtual void dkfdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dkfdc
            ) const = 0;

            //- Concentration derivative of reverse rate
            virtual void dkrdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalarField& dkfdc,
                const scalar kr,
                scalarField& dkrdc
            ) const = 0;

            //- Derivative of the net reaction rate for each species involved
            //  w.r.t. the concentration and temperature
            void ddNdtByVdcTp
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dNdtByV,
                scalarSquareMatrix& ddNdtByVdcTp,
                const bool reduced,
                const List<label>& c2s,
                const label csi0,
                const label Tsi,
                scalarField& cTpWork0,
                scalarField& cTpWork1
            ) const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Reaction<ThermoType>&) = delete;


    // Ostream Operator

        friend Ostream& operator<< <ThermoType>
        (
            Ostream&,
            const Reaction<ThermoType>&
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
