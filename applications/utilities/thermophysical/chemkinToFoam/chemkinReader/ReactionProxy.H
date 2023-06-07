/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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
    Foam::ReactionProxy

Description
    Proxy version of Reaction which provides dummy implementations of the
    abstract virtual functions.

    Used for read-construction and format conversion.

SourceFiles
    ReactionProxy.C

\*---------------------------------------------------------------------------*/

#ifndef ReactionProxy_H
#define ReactionProxy_H

#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ReactionProxy Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class ReactionProxy
:
    public Reaction<ThermoType>
{

public:

    // Constructors

        //- Construct from components
        ReactionProxy
        (
            const speciesTable& species,
            const PtrList<ThermoType>& speciesThermo,
            const List<specieCoeffs>& lhs,
            const List<specieCoeffs>& rhs
        );

        //- Construct and return a clone
        virtual autoPtr<Reaction<ThermoType>> clone() const;

        //- Construct and return a clone with new speciesTable
        virtual autoPtr<Reaction<ThermoType>> clone
        (
            const speciesTable& species
        ) const;


    //- Destructor
    virtual ~ReactionProxy()
    {}


    // Member Functions

        // Hooks

            //- Pre-evaluation hook
            virtual void preEvaluate() const;

            //- Post-evaluation hook
            virtual void postEvaluate() const;


        // Reaction rate coefficients

            //- Forward rate constant
            virtual scalar kf
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Reverse rate constant from the given forward rate constant
            virtual scalar kr
            (
                const scalar kfwd,
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Reverse rate constant
            virtual scalar kr
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;


        // Jacobian coefficients

            //- Temperature derivative of forward rate
            virtual scalar dkfdT
            (
                const scalar p,
                const scalar T,
                const scalarField& c,
                const label li
            ) const;

            //- Temperature derivative of reverse rate
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

            //- Concentration derivative of reverse rate
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactionProxy.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
