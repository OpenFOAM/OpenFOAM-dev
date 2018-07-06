/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::laminarModels::Maxwell

Description
    Maxwell model for viscoelasticity using the upper-convected time
    derivative of the stress tensor.
    See http://en.wikipedia.org/wiki/Upper-convected_Maxwell_model

    The model includes an additional viscosity (nu) from the transport
    model from which it is instantiated, which makes it equivalent to
    the Oldroyd-B model for the case of an incompressible transport
    model (where nu is non-zero).
    See https://en.wikipedia.org/wiki/Oldroyd-B_model

    Reference:
    \verbatim
        Amoreira, L. J., & Oliveira, P. J. (2010).
        Comparison of different formulations for the numerical calculation
        of unsteady incompressible viscoelastic fluid flow.
        Adv. Appl. Math. Mech, 4, 483-502.
    \endverbatim

SourceFiles
    Maxwell.C

\*---------------------------------------------------------------------------*/

#ifndef Maxwell_H
#define Maxwell_H

#include "laminarModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{

/*---------------------------------------------------------------------------*\
                           Class Maxwell Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class Maxwell
:
    public laminarModel<BasicTurbulenceModel>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar nuM_;
            dimensionedScalar lambda_;


        // Fields

            volSymmTensorField sigma_;


    // Protected Member Functions

        //- Return the turbulence viscosity
        tmp<volScalarField> nu0() const
        {
            return this->nu() + nuM_;
        }


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("Maxwell");


    // Constructors

        //- Construct from components
        Maxwell
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~Maxwell()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct eddy-Viscosity and
        //  related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Maxwell.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
