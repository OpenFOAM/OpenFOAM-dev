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
    Foam::blendingMethods::segregated

Description
    Blending method for segregated configurations. E.g., a churning flow in
    which it is never appropriate to consider one phase continuous. Applies
    only the segregated (i.e., phaseA_segregatedWith_phaseB) model.

    Example usage, in constant/phaseProperties:
    \verbatim
    blending
    {
        default
        {
            type            segregated;
        }
    }
    \endverbatim

SourceFiles
    segregated_blendingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef segregated_blendingMethod_H
#define segregated_blendingMethod_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                         Class segregated Declaration
\*---------------------------------------------------------------------------*/

class segregated
:
    public blendingMethod
{
protected:

    // Protected Member Functions

        //- Return the coefficient for models for which the indexed phase is
        //  considered continuous
        virtual tmp<volScalarField> fContinuous
        (
            const UPtrList<const volScalarField>& alphas,
            const label index
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
    TypeName("segregated");


    // Constructors

        //- Construct from a dictionary and an interface
        segregated
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        );


    //- Destructor
    ~segregated();


    // Member Functions

        //- Return whether or not a phase can be considered continuous
        virtual bool canBeContinuous(const label index) const;

        //- Return whether or not this interface can segregate
        virtual bool canSegregate() const;

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
