/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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
    Foam::phasePair

Description

SourceFiles
    phasePair.C

\*---------------------------------------------------------------------------*/

#ifndef phasePair_H
#define phasePair_H

#include "phaseModel.H"
#include "phasePairKey.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class phasePair Declaration
\*---------------------------------------------------------------------------*/

class phasePair
:
    public phasePairKey
{
private:

    // Private Data

        //- Phase 1
        const phaseModel& phase1_;

        //- Phase 2
        const phaseModel& phase2_;

        //- Gravitational acceleration
        const uniformDimensionedVectorField& g_;


public:

    // Constructors

        //- Construct from two phases and gravity
        phasePair
        (
            const phaseModel& phase1,
            const phaseModel& phase2,
            const bool ordered = false
        );


    //- Destructor
    virtual ~phasePair();


    // Member Functions

        //- Dispersed phase
        virtual const phaseModel& dispersed() const;

        //- Continuous phase
        virtual const phaseModel& continuous() const;

        //- Pair name
        virtual word name() const;

        //- Other pair name
        virtual word otherName() const;

        //- Average density
        tmp<volScalarField> rho() const;

        //- Relative velocity magnitude
        tmp<volScalarField> magUr() const;

        //- Relative velocity
        tmp<volVectorField> Ur() const;

        //- Reynolds number
        tmp<volScalarField> Re() const;

        //- Prandtl number
        tmp<volScalarField> Pr() const;

        //- Eotvos number
        tmp<volScalarField> Eo() const;

        //- Eotvos number for given diameter
        tmp<volScalarField> Eo(const volScalarField& d) const;

        //- Surface tension coefficient
        tmp<volScalarField> sigma() const;

        //- Morton Number
        tmp<volScalarField> Mo() const;

        //- Takahashi Number
        tmp<volScalarField> Ta() const;


        // Access

            //- Return phase 1
            inline const phaseModel& phase1() const;

            //- Return phase 2
            inline const phaseModel& phase2() const;

            //- Return true if this phasePair contains the given phase
            inline bool contains(const phaseModel& phase) const;

            //- Return the other phase relative to the given phase
            //  Generates a FatalError if this phasePair does not contain
            //  the given phase
            inline const phaseModel& otherPhase(const phaseModel& phase) const;

            //- Return the index of the given phase. Generates a FatalError if
            //  this phasePair does not contain the given phase
            inline label index(const phaseModel& phase) const;

            //- Return gravitation acceleration
            inline const uniformDimensionedVectorField& g() const;

            //- Return the mesh
            inline const fvMesh& mesh() const;


        //- STL const_iterator
        class const_iterator
        {
            // Private Data

                //- Reference to the pair for which this is an iterator
                const phasePair& pair_;

                //- Current index
                label index_;

                //- Construct an iterator with the given index
                inline const_iterator(const phasePair&, const label index);

        public:

            friend class phasePair;

            // Constructors

                //- Construct from pair, moving to its 'begin' position
                inline explicit const_iterator(const phasePair&);


            // Access

                //- Return the current index
                inline label index() const;


            // Member Operators

                inline bool operator==(const const_iterator&) const;

                inline bool operator!=(const const_iterator&) const;

                inline const phaseModel& operator*() const;
                inline const phaseModel& operator()() const;

                inline const phaseModel& otherPhase() const;

                inline const_iterator& operator++();
                inline const_iterator operator++(int);
        };


        //- const_iterator set to the beginning of the pair
        inline const_iterator cbegin() const;

        //- const_iterator set to beyond the end of the pair
        inline const_iterator cend() const;

        //- const_iterator set to the beginning of the pair
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end of the pair
        inline const_iterator end() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phasePairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
