/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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
    Foam::RASModel

Description
    Templated abstract base class for RAS turbulence models

    with support for generalised Newtonian viscosity models including
    strain-rate dependency.

SourceFiles
    RASModel.C

\*---------------------------------------------------------------------------*/

#ifndef RASModel_H
#define RASModel_H

#include "momentumTransportModel.H"
#include "generalisedNewtonianViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class RASModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class RASModel
:
    public BasicMomentumTransportModel
{

protected:

    // Protected data

        //- RAS coefficients dictionary
        dictionary RASDict_;

        //- Turbulence on/off flag
        Switch turbulence_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- Lower limit of k
        dimensionedScalar kMin_;

        //- Lower limit of epsilon
        dimensionedScalar epsilonMin_;

        //- Lower limit for omega
        dimensionedScalar omegaMin_;

        //- Run-time selectable generalised Newtonian viscosity model
        autoPtr<laminarModels::generalisedNewtonianViscosityModel>
            viscosityModel_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs(const word& type);


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("RAS");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RASModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const viscosity& viscosity
            ),
            (alpha, rho, U, alphaRhoPhi, phi, viscosity)
        );


    // Constructors

        //- Construct from components
        RASModel
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );

        //- Disallow default bitwise copy construction
        RASModel(const RASModel&) = delete;


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<RASModel> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~RASModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the lower allowable limit for k (default: small)
        const dimensionedScalar& kMin() const
        {
            return kMin_;
        }

        //- Return the lower allowable limit for epsilon (default: small)
        const dimensionedScalar& epsilonMin() const
        {
            return epsilonMin_;
        }

        //- Return the lower allowable limit for omega (default: small)
        const dimensionedScalar& omegaMin() const
        {
            return omegaMin_;
        }

        //- Allow kMin to be changed
        dimensionedScalar& kMin()
        {
            return kMin_;
        }

        //- Allow epsilonMin to be changed
        dimensionedScalar& epsilonMin()
        {
            return epsilonMin_;
        }

        //- Allow omegaMin to be changed
        dimensionedScalar& omegaMin()
        {
            return omegaMin_;
        }

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return viscosityModel_->nu();
        }

        //- Return the laminar viscosity on patchi
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return viscosityModel_->nu(patchi);
        }

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            return volScalarField::New
            (
                this->groupName("nuEff"),
                this->nut() + this->nu()
            );
        }

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const
        {
            return this->nut(patchi) + this->nu(patchi);
        }

        //- Predict the turbulence transport coefficients if possible
        //  without solving turbulence transport model equations
        virtual void predict()
        {}

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const RASModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "RASModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
