/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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
    Foam::blendingMethods::continuous

Description
    Blending method for the case in which the continuous phase is always the
    same. E.g., for an air-particles flow, the air phase is always continuous.
    Applies only the dispersed (i.e., phaseA_dispersedIn_phaseB, where phaseB
    is the identified continuous phase) model.

    Example usage, in constant/phaseProperties:
    \verbatim
    blending
    {
        default
        {
            type            continuous;
            phase           air;
        }
    }
    \endverbatim

SourceFiles
    continuous.C

\*---------------------------------------------------------------------------*/

#ifndef continuous_H
#define continuous_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                         Class continuous Declaration
\*---------------------------------------------------------------------------*/

class continuous
:
    public blendingMethod
{
    // Private Data

        //- The continuous phase
        const phaseModel& phase_;

        //- Apply displaced models? The displaced model coefficient will be
        //  calculated as one-minus the continuous coefficient calculated for
        //  the two-phase sub-system.
        const bool displaced_;

        //- The index of the continuous phase in the interface, or -1 if the
        //  interface does not contain the continuous phase
        const label index_;

        //- The index of the non-continuous phase in the interface, or -1 if
        //  the interface does not contain the continuous phase
        const label otherIndex_;


protected:

    // Protected Member Functions

        //- Return the coefficient for models for which the indexed phase is
        //  considered continuous
        virtual tmp<volScalarField> fContinuous
        (
            const UPtrList<const volScalarField>& alphas,
            const label index
        ) const;

        //- Return the coefficient for segregated models
        virtual tmp<volScalarField> fSegregated
        (
            const UPtrList<const volScalarField>& alphas
        ) const;

        //- Return the coefficient for when the interface is displaced by a
        //  third phase
        virtual tmp<volScalarField> fDisplaced
        (
            const UPtrList<const volScalarField>& alphas,
            const label displacingPhasei
        ) const;


public:

    //- Runtime type information
    TypeName("continuous");


    // Constructors

        //- Construct from a dictionary and an interface
        continuous
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        );


    //- Destructor
    ~continuous();


    // Member Functions

        //- Return whether or not a phase can be considered continuous
        virtual bool canBeContinuous(const label index) const;

        //- Return whether or not this interface can segregate
        virtual bool canSegregate() const;

        //- Return whether or not this method partitions unity across the
        //  entire blending space
        virtual bool isComplete() const;

        //- Return whether or not the interface is displaced by a third phase
        virtual bool isDisplacedBy(const label displacingPhasei) const;

        //- Return whether this is a function of the volume fractions only
        virtual bool functionOfAlphas() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
