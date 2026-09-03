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
    Foam::blendingMethods::linearField

Description
    Blending method which combines two sub-blending methods based on a linear
    function of a user-specified field. Two values are specified which defines
    the range of field values over which the blending is transitioned between
    the two sub-blending approaches.

    Example usage, in constant/phaseProperties:
    \verbatim
    blending
    {
        default
        {
            type            linearField;
            field           myField;
            values          (0.4 0.6);
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
    linearField.C

\*---------------------------------------------------------------------------*/

#ifndef linearField_H
#define linearField_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                           Class linearField Declaration
\*---------------------------------------------------------------------------*/

class linearField
:
    public blendingMethod
{
    // Private Data

        //- The name of the field which defines the blending regions
        const word fieldName_;

        //- Coefficients dictionary held temporarily until the field value can
        //  be constructed
        mutable autoPtr<dictionary> dictPtr_;

        //- The values of the field between which the blending switches
        mutable autoPtr<Pair<dimensionedScalar>> fieldValuesPtr_;

        //- The blending when the field is below the specified value
        autoPtr<blendingMethod> belowBlending_;

        //- The blending when the field is above the specified value
        autoPtr<blendingMethod> aboveBlending_;


    // Private Member Functions

        //- Return the value of the field at which the blending switches
        const Pair<dimensionedScalar>& fieldValues() const;


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
    TypeName("linearField");


    // Constructors

        //- Construct from a dictionary and an interface
        linearField
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        );


    //- Destructor
    ~linearField();


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
