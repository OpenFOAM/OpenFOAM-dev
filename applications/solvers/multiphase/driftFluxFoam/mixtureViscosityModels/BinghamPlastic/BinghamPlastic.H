/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::mixtureViscosityModels::BinghamPlastic

Description
     Viscosity correction model for Bingham plastics.

    The strain-rate used is defined as sqrt(2.0)*mag(symm(grad(U)))

SourceFiles
    BinghamPlastic.C

\*---------------------------------------------------------------------------*/

#ifndef BinghamPlastic_H
#define BinghamPlastic_H

#include "plastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class BinghamPlastic Declaration
\*---------------------------------------------------------------------------*/

class BinghamPlastic
:
    public plastic
{
protected:

    // Protected data

        //- Yield stress coefficient
        dimensionedScalar yieldStressCoeff_;

        //- Yield stress exponent
        dimensionedScalar yieldStressExponent_;

        //- Yield stress offset
        dimensionedScalar yieldStressOffset_;

        //- Velocity
        const volVectorField& U_;


public:

    //- Runtime type information
    TypeName("BinghamPlastic");


    // Constructors

        //- Construct from components
        BinghamPlastic
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~BinghamPlastic()
    {}


    // Member Functions

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        tmp<volScalarField> mu(const volScalarField& muc) const;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
