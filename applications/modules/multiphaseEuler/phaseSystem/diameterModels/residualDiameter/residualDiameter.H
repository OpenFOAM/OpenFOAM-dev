/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::diameterModels::residualDiameter

Description
    A diameter model which switches form constant diameter to constant residual
    diameter when the volume fraction of the phase is below residualAlpha.

    A typical use case for this is to set a small residual diameter to prevent
    excessive phase velocities in regions where the phase volume fraction is
    small.

SourceFiles
    residualDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef residualDiameter_H
#define residualDiameter_H

#include "sphericalDiameter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                      Class residualDiameter Declaration
\*---------------------------------------------------------------------------*/

class residualDiameter
:
    public spherical
{
    // Private data

        //- The diameter used when the volume fraction is above residualAlpha
        dimensionedScalar d_;

        //- The diameter used when the volume fraction is below residualAlpha
        dimensionedScalar dResidual_;


public:

    //- Runtime type information
    TypeName("residualDiameter");


    // Constructors

        residualDiameter
        (
            const dictionary& diameterProperties,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~residualDiameter();


    // Member Functions

        //- Get the diameter field
        virtual tmp<volScalarField> d() const;

        //- Read diameterProperties dictionary
        virtual bool read(const dictionary& diameterProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
