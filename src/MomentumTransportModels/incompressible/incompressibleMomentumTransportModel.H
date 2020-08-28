/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::incompressibleMomentumTransportModel

Description
    Abstract base class for turbulence models (RAS, LES and laminar).

SourceFiles
    incompressibleMomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleMomentumTransportModel_H
#define incompressibleMomentumTransportModel_H

#include "momentumTransportModel.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
            Class incompressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class incompressibleMomentumTransportModel
:
    public momentumTransportModel
{
protected:

    // Protected data

        geometricOneField rho_;


public:

    //- Runtime type information
    TypeName(momentumTransportModel::typeName_());


    // Constructors

        //- Construct from components
        incompressibleMomentumTransportModel
        (
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        incompressibleMomentumTransportModel
        (
            const incompressibleMomentumTransportModel&
        ) = delete;


    //- Destructor
    virtual ~incompressibleMomentumTransportModel()
    {}


    // Member Functions

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const;

        //- Return the laminar dynamic viscosity on patch
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Return the turbulence dynamic viscosity
        virtual tmp<volScalarField> mut() const;

        //- Return the turbulence dynamic viscosity on patch
        virtual tmp<scalarField> mut(const label patchi) const;

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const;

        //- Return the effective dynamic viscosity on patch
        virtual tmp<scalarField> muEff(const label patchi) const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devSigma() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevSigma(volVectorField& U) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const incompressibleMomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
