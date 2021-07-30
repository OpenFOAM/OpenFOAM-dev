/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
    Foam::laminarModels::generalisedNewtonianViscosityModels::
    strainRateViscosityModel

Description
    An abstract base class for strain-rate dependent generalised Newtonian
    viscosity models

SourceFiles
    strainRateViscosityModel.C

\*---------------------------------------------------------------------------*/

#ifndef strainRateViscosityModel_H
#define strainRateViscosityModel_H

#include "generalisedNewtonianViscosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                Class strainRateViscosityModel Declaration
\*---------------------------------------------------------------------------*/

class strainRateViscosityModel
:
    public generalisedNewtonianViscosityModel
{
    // Private data

        //- The strain-rate-dependent generalised Newtonian viscosity field
        volScalarField nu_;


protected:

    // Protected Member Functions

        //- Return the strain-rate obtained from the velocity field
        virtual tmp<volScalarField> strainRate() const;

        //- Return the generalised Newtonian viscosity
        virtual tmp<volScalarField> nu
        (
            const volScalarField& nu0,
            const volScalarField& strainRate
        ) const = 0;


public:

    //- Runtime type information
    TypeName("strainRateViscosityModel");


    // Constructors

        //- Construct from components
        strainRateViscosityModel
        (
            const dictionary& viscosityProperties,
            const viscosity& viscosity,
            const volVectorField& U
        );

        //- Disallow default bitwise copy construction
        strainRateViscosityModel
        (
            const strainRateViscosityModel&
        ) = delete;


    //- Destructor
    virtual ~strainRateViscosityModel()
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
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const strainRateViscosityModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace generalisedNewtonianViscosityModels
} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
