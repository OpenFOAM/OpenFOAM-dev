/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::RASModels::kOmega2006

Description
    Standard (2006) high Reynolds-number k-omega turbulence model for
    incompressible and compressible flows.

    References:
    \verbatim
        Wilcox, D. C. (2006).
        Turbulence modeling for CFD, 3rd edition
        La Canada, CA: DCW industries, Inc.

        Wilcox, D. C. (2008).
        Formulation of the kw turbulence model revisited.
        AIAA journal, 46(11), 2823-2838.
    \endverbatim

    The default model coefficients are
    \verbatim
        kOmega2006Coeffs
        {
            Cmu         0.09;
            beta0       0.0708;
            gamma       0.52;
            Clim        0.875;
            alphak      0.6;
            alphaOmega  0.5;
        }
    \endverbatim

SourceFiles
    kOmega2006.C

\*---------------------------------------------------------------------------*/

#ifndef kOmega2006_H
#define kOmega2006_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmega2006 Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmega2006
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
    // Private Member Functions

        using IOobject::modelName;

        void correctNut(const volTensorField& gradU);

        tmp<volScalarField::Internal> beta(const volTensorField& gradU) const;

        tmp<volScalarField::Internal> CDkOmega() const;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar betaStar_;
            dimensionedScalar beta0_;
            dimensionedScalar gamma_;
            dimensionedScalar Clim_;
            dimensionedScalar sigmaDo_;
            dimensionedScalar alphaK_;
            dimensionedScalar alphaOmega_;


        // Fields

            volScalarField k_;
            volScalarField omega_;


    // Protected Member Functions

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> omegaSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmega2006");


    // Constructors

        //- Construct from components
        kOmega2006
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmega2006()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                alphaK_*this->nut_ + this->nu()
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return volScalarField::New
            (
                "DomegaEff",
                alphaOmega_*this->nut_ + this->nu()
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return volScalarField::New
            (
                "epsilon",
                betaStar_*k_*omega_,
                omega_.boundaryField().types()
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmega2006.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
