/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::viscosityModels::CrossPowerLaw

Description
    An incompressible Cross-Power law non-Newtonian viscosity model.

SourceFiles
    CrossPowerLaw.C

\*---------------------------------------------------------------------------*/

#ifndef CrossPowerLaw_H
#define CrossPowerLaw_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                        Class CrossPowerLaw Declaration
\*---------------------------------------------------------------------------*/

class CrossPowerLaw
:
    public viscosityModel
{
    // Private Data

        dictionary CrossPowerLawCoeffs_;

        dimensionedScalar nu0_;
        dimensionedScalar nuInf_;
        dimensionedScalar m_;
        dimensionedScalar n_;

        volScalarField nu_;

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("CrossPowerLaw");


    // Constructors

        //- Construct from components
        CrossPowerLaw
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~CrossPowerLaw()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
