/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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
    Foam::mixtureViscosityModels::plastic

Description
     Viscosity correction model for a generic power-law plastic.

SourceFiles
    plastic.C

\*---------------------------------------------------------------------------*/

#ifndef plastic_H
#define plastic_H

#include "mixtureViscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleDriftFluxMixture;

namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class plastic Declaration
\*---------------------------------------------------------------------------*/

class plastic
:
    public mixtureViscosityModel
{
protected:

    // Protected data

        //- Plastic viscosity coefficient
        dimensionedScalar plasticViscosityCoeff_;

        //- Plastic viscosity exponent
        dimensionedScalar plasticViscosityExponent_;

        //- Maximum viscosity
        dimensionedScalar muMax_;


public:

    //- Runtime type information
    TypeName("plastic");


    // Constructors

        //- Construct from mixture
        plastic(const incompressibleDriftFluxMixture& mixture);


    //- Destructor
    virtual ~plastic()
    {}


    // Member Functions

        const dictionary& coeffDict() const
        {
            return optionalSubDict(typeName + "Coeffs");
        }

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        virtual tmp<volScalarField> mu
        (
            const volScalarField& muc,
            const volVectorField& U
        ) const;

        //- Read phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
