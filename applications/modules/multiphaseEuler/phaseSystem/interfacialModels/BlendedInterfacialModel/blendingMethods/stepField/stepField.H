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
    Foam::blendingMethods::stepField

Description
    Blending method which switches between two sub-blending methods based on
    the value of a user-specified field. A values is specified which defines
    the field value at which the blending switches between the two sub-blending
    approaches.

    Example usage, in constant/phaseProperties:
    \verbatim
    blending
    {
        default
        {
            type            stepField;
            field           myField;
            value           0.5;
            below
            {
                type            linear;
                displaced       no;
                minFullyContinuousAlpha.air 0.8;
                minPartlyContinuousAlpha.air 0.2;
                minFullyContinuousAlpha.water 0.8;
                minPartlyContinuousAlpha.water 0.2;
                minFullyContinuousAlpha.particles none;
                minPartlyContinuousAlpha.particles none;
                minFullyContinuousAlpha.bed none;
                minPartlyContinuousAlpha.bed none;
            }
            above
            {
                type            segregated;
            }
        }
    }
    \endverbatim

SourceFiles
    stepField.C

\*---------------------------------------------------------------------------*/

#ifndef stepField_H
#define stepField_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                           Class stepField Declaration
\*---------------------------------------------------------------------------*/

class stepField
:
    public blendingMethod
{
    // Private Data

        //- The name of the field which defines the blending regions
        const word fieldName_;

        //- Coefficients dictionary held temporarily until the field value can
        //  be constructed
        mutable autoPtr<dictionary> dictPtr_;

        //- The value of the field at which the blending switches
        mutable autoPtr<dimensionedScalar> fieldValuePtr_;

        //- The blending when the field is below the specified value
        autoPtr<blendingMethod> belowBlending_;

        //- The blending when the field is above the specified value
        autoPtr<blendingMethod> aboveBlending_;


    // Private Member Functions

        //- Return the value of the field at which the blending switches
        const dimensionedScalar& fieldValue() const;


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
    TypeName("stepField");


    // Constructors

        //- Construct from a dictionary and an interface
        stepField
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        );


    //- Destructor
    ~stepField();


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
