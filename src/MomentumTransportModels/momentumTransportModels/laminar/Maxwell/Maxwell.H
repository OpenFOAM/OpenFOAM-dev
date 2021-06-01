/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Generalised Maxwell model for viscoelasticity using the upper-convected time
    derivative of the stress tensor with support for multiple modes.

    See http://en.wikipedia.org/wiki/Upper-convected_Maxwell_model
        http://en.wikipedia.org/wiki/Generalised_Maxwell_model

    The model includes an additional viscosity (nu) from the transport
    model from which it is instantiated, which makes it equivalent to
    the Oldroyd-B model for the case of an incompressible transport
    model (where nu is non-zero).
    See https://en.wikipedia.org/wiki/Oldroyd-B_model

    Reference:
    \verbatim
        Wiechert, E. (1889). Ueber elastische Nachwirkung.
        (Doctoral dissertation, Hartungsche buchdr.).

        Wiechert, E. (1893).
        Gesetze der elastischen Nachwirkung für constante Temperatur.
        Annalen der Physik, 286(11), 546-570.

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

template<class BasicMomentumTransportModel>
class Maxwell
:
    public laminarModel<BasicMomentumTransportModel>
{

protected:

    // Protected data

        // Model coefficients

            PtrList<dictionary> modeCoefficients_;

            label nModes_;

            dimensionedScalar nuM_;

            PtrList<dimensionedScalar> lambdas_;


        // Fields

            //- Single or mode sum viscoelastic stress
            volSymmTensorField sigma_;

            //- Mode viscoelastic stresses
            PtrList<volSymmTensorField> sigmas_;


    // Protected Member Functions

        PtrList<dimensionedScalar> readModeCoefficients
        (
            const word& name,
            const dimensionSet& dims
        ) const;

        //- Return the turbulence viscosity
        tmp<volScalarField> nu0() const
        {
            return this->nu() + nuM_;
        }

        virtual tmp<fvSymmTensorMatrix> sigmaSource
        (
            const label modei,
            volSymmTensorField& sigma
        ) const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;


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
            const word& type = typeName
        );


    //- Destructor
    virtual ~Maxwell()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the effective viscosity, i.e. the laminar viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective viscosity on patch
        virtual tmp<scalarField> nuEff(const label patchi) const;

        //- Return the stress tensor [m^2/s^2]
        virtual tmp<volSymmTensorField> sigma() const;

        //- Return the effective stress tensor
        virtual tmp<volSymmTensorField> devTau() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau
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
