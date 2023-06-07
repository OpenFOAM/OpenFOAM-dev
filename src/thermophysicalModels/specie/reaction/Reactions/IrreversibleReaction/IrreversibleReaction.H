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
    Foam::IrreversibleReaction

Description
    Extension of Reaction to handle irreversible reactions

SourceFiles
    IrreversibleReaction.C

\*---------------------------------------------------------------------------*/

#ifndef IrreversibleReaction_H
#define IrreversibleReaction_H

#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class IrreversibleReaction Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType, class ReactionRate>
class IrreversibleReaction
:
    public Reaction<ThermoType>
{
    // Private Data

        ReactionRate k_;


public:

    //- Runtime type information
    TypeName("irreversible");


    // Constructors

        //- Construct from components
        IrreversibleReaction
        (
            const Reaction<ThermoType>& reaction,
            const ReactionRate& reactionRate
        );

        //- Construct as copy given new speciesTable
        IrreversibleReaction
        (
            const IrreversibleReaction<ThermoType, ReactionRate>&,
            const speciesTable& species
        );

        //- Construct from dictionary
        IrreversibleReaction
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const dictionary& dict
        );

        //- Construct from objectRegistry and dictionary
        IrreversibleReaction
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const objectRegistry& ob,
            const dictionary& dict
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ThermoType>> clone() const
        {
            return autoPtr<Reaction<ThermoType>>
            (
                new IrreversibleReaction<ThermoType, ReactionRate>
                (
                    *this
                )
            );
        }

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ThermoType>> clone
        (
            const speciesTable& species
        ) const
        {
            return autoPtr<Reaction<ThermoType>>
            (
                new IrreversibleReaction<ThermoType, ReactionRate>
                (
                    *this,
                    species
                )
            );
        }


    //- Destructor
    virtual ~IrreversibleReaction()
    {}


    // Member Functions

        // Hooks

            //- Pre-evaluation hook
            virtual void preEvaluate() const;

            //- Post-evaluation hook
            virtual void postEvaluate() const;


        // IrreversibleReaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Reverse rate constant from the given forward rate constant
            //  Returns 0
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Reverse rate constant
            //  Returns 0
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;


        // IrreversibleReaction Jacobian functions

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Temperature derivative of reverse rate. Returns zero.
            virtual scalar dkrdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalar dkfdT,
                const scalar kr
            ) const;

            //- Does this reaction have concentration-dependent rate constants?
            virtual bool hasDkdc() const;

            //- Concentration derivative of forward rate
            void dkfdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                scalarField& dkfdc
            ) const;

            //- Concentration derivative of reverse rate. Returns zero.
            void dkrdc
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li,
                const scalarField& dkfdc,
                const scalar kr,
                scalarField& dkrdc
            ) const;


        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=
        (
            const IrreversibleReaction<ThermoType, ReactionRate>&
        ) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IrreversibleReaction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
