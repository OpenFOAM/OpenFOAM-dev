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
    Foam::blendingMethod

Description
    Abstract base class for functions that are used to combine interfacial
    sub-models according to the volume fractions of the phases that they apply
    to

SourceFiles
    blendingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef blendingMethod_H
#define blendingMethod_H

#include "phaseInterface.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Struct blendingParameter Declaration
\*---------------------------------------------------------------------------*/

struct blendingParameter
{
    //- Was this parameter specified?
    bool specified;

    //- The parameter value. Could be templated if needed.
    scalar value;
};


/*---------------------------------------------------------------------------*\
                       Class blendingMethod Declaration
\*---------------------------------------------------------------------------*/

class blendingMethod
{
private:

    // Private member functions

        //- Read a parameter and check it lies within specified bounds
        static blendingParameter readParameter
        (
            const word& name,
            const dictionary& dict,
            const Pair<scalar>& bounds,
            const bool allowNone,
            const scalar noneValue
        );


protected:

    // Protected Data

        //- Interface
        const phaseInterface interface_;


    // Protected Static Member Functions

        //- Read a parameter and check it lies within specified bounds
        static blendingParameter readParameter
        (
            const word& name,
            const dictionary& dict,
            const Pair<scalar>& bounds
        );

        //- Read a parameter and check it lies within specified bounds
        static blendingParameter readParameter
        (
            const word& name,
            const dictionary& dict,
            const Pair<scalar>& bounds,
            const scalar noneValue
        );

        //- Read a parameter and check it lies within specified bounds
        static blendingParameter readParameterIfPresent
        (
            const word& name,
            const dictionary& dict,
            const Pair<scalar>& bounds,
            const scalar noneValue
        );

        //- Read a parameter for each phase in the interface
        static Pair<blendingParameter> readParameters
        (
            const word& name,
            const dictionary& dict,
            const phaseInterface& interface,
            const Pair<scalar>& bounds
        );

        //- Read a parameter for each phase in the interface
        static Pair<blendingParameter> readParameters
        (
            const word& name,
            const dictionary& dict,
            const phaseInterface& interface,
            const Pair<scalar>& bounds,
            const scalar noneValue
        );


    // Protected Member Functions

        //- Return a constant field with the given value
        tmp<volScalarField> constant
        (
            const UPtrList<const volScalarField>& alphas,
            const scalar k
        ) const;

        //- Get the volume fraction given the index in the interface. If the
        //  index is 0 or 1 this returns the volume fraction of the indexed
        //  phase, and if the index is -1 then this returns the volume fraction
        //  of the two-phase sub-system.
        tmp<volScalarField> alpha
        (
            const UPtrList<const volScalarField>& alphas,
            const label index,
            const bool protect
        ) const;

        //- Get a blending parameter. Index has the same meaning as for
        //  blendingMethod::alpha. If the index is -1 then the parameter
        //  returned will be a value constructed for the sub-system by
        //  volume-fraction-weighting the parameter values for both phases.
        tmp<volScalarField> parameter
        (
            const UPtrList<const volScalarField>& alphas,
            const label index,
            const Pair<blendingParameter>& parameters
        ) const;

        //- Return the coordinate of the blending function. Index has the same
        //  meaning as for blendingMethod::alpha.
        tmp<volScalarField> x
        (
            const UPtrList<const volScalarField>& alphas,
            const label index
        ) const;

        //- Return the coefficient for models for which the indexed phase is
        //  considered continuous. Index has the same meaning as for
        //  blendingMethod::alpha.
        virtual tmp<volScalarField> fContinuous
        (
            const UPtrList<const volScalarField>& alphas,
            const label index
        ) const = 0;

        //- As above, but for a different blending method
        static inline tmp<volScalarField> fContinuous
        (
            const blendingMethod& blending,
            const UPtrList<const volScalarField>& alphas,
            const label index
        )
        {
            return blending.fContinuous(alphas, index);
        }

        //- Return the coefficient for when the interface is displaced by a
        //  third phase
        virtual tmp<volScalarField> fDisplaced
        (
            const UPtrList<const volScalarField>& alphas,
            const label displacingPhasei
        ) const = 0;

        //- As above, but for a different blending method
        static inline tmp<volScalarField> fDisplaced
        (
            const blendingMethod& blending,
            const UPtrList<const volScalarField>& alphas,
            const label displacingPhasei
        )
        {
            return blending.fDisplaced(alphas, displacingPhasei);
        }

        //- As blendingMethod::fContinuous, but filters out and returns default
        //  values if the indexed phase cannot become continuous, then calls
        //  the virtual fContinuous.
        tmp<volScalarField> fContinuousFiltered
        (
            const UPtrList<const volScalarField>& alphas,
            const label index
        ) const;


public:

    //- Runtime type information
    TypeName("blendingMethod");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        blendingMethod,
        dictionary,
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced
        ),
        (dict, interface, allowDisplaced)
    );


    // Constructors

        //- Construct for an interface
        blendingMethod(const phaseInterface& interface);


    // Selector

        static autoPtr<blendingMethod> New
        (
            const word& modelTypeName,
            const dictionary& dict,
            const phaseInterface& interface,
            const bool allowDisplaced = true
        );


    //- Destructor
    virtual ~blendingMethod();


    // Member Functions

        //- Return whether or not a phase can be considered continuous
        virtual bool canBeContinuous(const label index) const = 0;

        //- Return whether or not this interface can segregate
        virtual bool canSegregate() const = 0;

        //- Return whether or not the interface is displaced by a third phase
        virtual bool isDisplacedBy(const label displacingPhasei) const = 0;

        //- Return whether this is a function of the volume fractions only
        virtual bool functionOfAlphas() const = 0;

        //- Return the coefficient for models in which phase 1 is dispersed in
        //  phase 2
        tmp<volScalarField> f1DispersedIn2
        (
            const UPtrList<const volScalarField>& alphas
        ) const;

        //- Return the coefficient for models in which phase 2 is dispersed in
        //  phase 1
        tmp<volScalarField> f2DispersedIn1
        (
            const UPtrList<const volScalarField>& alphas
        ) const;

        //- Return the coefficient for models in which phases 1 and 2 are
        //  displaced by a third phase
        tmp<volScalarField> f12DisplacedBy3
        (
            const UPtrList<const volScalarField>& alphas,
            const label displacingPhasei
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
