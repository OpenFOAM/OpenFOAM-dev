/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

    This method requires the continuous phase to be specified.

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


protected:

    // Protected Member Functions

        //- Evaluate the blending function
        virtual tmp<volScalarField> fContinuous
        (
            const UPtrList<const volScalarField>& alphas,
            const label phaseSet,
            const label systemSet
        ) const;


public:

    //- Runtime type information
    TypeName("continuous");


    // Constructors

        //- Construct from a dictionary and an interface
        continuous
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    ~continuous();


    // Member Functions

        //- Return whether or not a phase can be considered continuous
        virtual bool canBeContinuous(const label index) const;

        //- Return whether or not this interface can segregate
        virtual bool canSegregate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
