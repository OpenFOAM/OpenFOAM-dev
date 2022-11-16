/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
    to.

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
                       Class blendingMethod Declaration
\*---------------------------------------------------------------------------*/

class blendingMethod
{
protected:

    // Protected Data

        //- Interface
        const phaseInterface interface_;


    // Protected Static Member Functions

        //- Read a parameter and check it lies within specified bounds
        static scalar readParameter
        (
            const word& name,
            const dictionary& dict,
            const Pair<scalar>& bounds,
            const bool allowNone
        );

        //- Read a parameter for each phase in the interface
        static Pair<scalar> readParameters
        (
            const word& name,
            const dictionary& dict,
            const phaseInterface& interface,
            const Pair<scalar>& bounds,
            const bool allowNone
        );

        //- Test if the given scalar is a valid parameter
        static bool isParameter(const scalar parameter);


    // Protected Member Functions

        //- Return a constant field with the given value
        tmp<volScalarField> constant
        (
            const UPtrList<const volScalarField>& alphas,
            const scalar k
        ) const;

        //- Get the volume fraction of the given set
        tmp<volScalarField> alpha
        (
            const UPtrList<const volScalarField>& alphas,
            const label set,
            const bool protect
        ) const;

        //- Get a blending parameter averaged for the given set
        tmp<volScalarField> parameter
        (
            const UPtrList<const volScalarField>& alphas,
            const label set,
            const Pair<scalar>& parameters
        ) const;

        //- Return the coordinate of the blending function
        tmp<volScalarField> x
        (
            const UPtrList<const volScalarField>& alphas,
            const label phaseSet,
            const label systemSet
        ) const;

        //- Evaluate the blending function for sets in which all phases can be
        //  continuous
        virtual tmp<volScalarField> fContinuous
        (
            const UPtrList<const volScalarField>& alphas,
            const label phaseSet,
            const label systemSet
        ) const = 0;

        //- Evaluate the blending function. Filters out phases that cannot
        //  become continuous from the sets, then calls fContinuous
        virtual tmp<volScalarField> f
        (
            const UPtrList<const volScalarField>& alphas,
            const label phaseSet,
            const label systemSet
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
            const phaseInterface& interface
        ),
        (dict, interface)
    );


    // Constructors

        //- Construct from a dictionary and an interface
        blendingMethod
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Selector

        static autoPtr<blendingMethod> New
        (
            const word& modelTypeName,
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~blendingMethod();


    // Member Functions

        //- Return whether or not a phase can be considered continuous
        virtual bool canBeContinuous(const label index) const = 0;

        //- Return whether or not this interface can segregate
        virtual bool canSegregate() const = 0;

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

        //- Return the coefficient for when the interface is displaced by a
        //  third phase
        tmp<volScalarField> fDisplaced
        (
            const UPtrList<const volScalarField>& alphas
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
