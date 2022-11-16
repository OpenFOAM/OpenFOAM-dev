/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::diameterModels::isothermal

Description
    Isothermal dispersed-phase particle diameter model.

SourceFiles
    isothermal.C

\*---------------------------------------------------------------------------*/

#ifndef isothermalDiameter_H
#define isothermalDiameter_H

#include "sphericalDiameter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                         Class isothermal Declaration
\*---------------------------------------------------------------------------*/

class isothermal
:
    public spherical
{
    // Private Data

        //- Reference diameter for the isothermal expansion
        dimensionedScalar d0_;

        //- Reference pressure for the isothermal expansion
        dimensionedScalar p0_;

        //- Diameter field
        volScalarField d_;


public:

    //- Runtime type information
    TypeName("isothermal");


    // Constructors

        //- Construct from dictionary and phase
        isothermal
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~isothermal();


    // Member Functions

        //- Get the diameter field
        virtual tmp<volScalarField> d() const;

        //- Correct the model
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
