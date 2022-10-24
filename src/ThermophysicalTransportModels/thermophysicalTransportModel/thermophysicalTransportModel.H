/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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
    Foam::thermophysicalTransportModel

Description
    Abstract base class for all fluid and solid thermophysical transport models

SourceFiles
    thermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef thermophysicalTransportModel_H
#define thermophysicalTransportModel_H

#include "IOdictionary.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class thermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

class thermophysicalTransportModel
:
    public IOdictionary
{

public:

    //- Runtime type information
    TypeName("thermophysicalTransport");


    // Constructors

        //- Construct from mesh and group
        thermophysicalTransportModel
        (
            const fvMesh& mesh,
            const word& group
        );

        //- Disallow default bitwise copy construction
        thermophysicalTransportModel
        (
            const thermophysicalTransportModel&
        ) = delete;


    //- Destructor
    virtual ~thermophysicalTransportModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read() = 0;

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const = 0;

        //- Effective thermal turbulent conductivity
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const = 0;

        //- Effective thermal turbulent conductivity
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const = 0;

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const = 0;

        //- Return the patch heat flux correction [W/m^2]
        //  For isotropic or patch-aligned thermal conductivity qCorr is null
        virtual tmp<scalarField> qCorr(const label patchi) const = 0;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const = 0;

        //- Solve the thermophysical transport model equations
        //  and correct the transport coefficients
        virtual void correct() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermophysicalTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
