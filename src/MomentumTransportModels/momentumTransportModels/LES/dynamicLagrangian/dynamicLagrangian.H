/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::LESModels::dynamicLagrangian

Description
    Dynamic SGS model with Lagrangian averaging

    Reference:
    \verbatim
        Meneveau, C., Lund, T. S., & Cabot, W. H. (1996).
        A Lagrangian dynamic subgrid-scale model of turbulence.
        Journal of Fluid Mechanics, 319, 353-385.
    \endverbatim

SourceFiles
    dynamicLagrangian.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicLagrangian_H
#define dynamicLagrangian_H

#include "LESeddyViscosity.H"
#include "simpleFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                      Class dynamicLagrangian Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class dynamicLagrangian
:
    public LESeddyViscosity<BasicMomentumTransportModel>
{
protected:

    // Protected data

        volScalarField flm_;
        volScalarField fmm_;

        dimensionedScalar theta_;

        simpleFilter simpleFilter_;
        autoPtr<LESfilter> filterPtr_;
        LESfilter& filter_;

        dimensionedScalar flm0_;
        dimensionedScalar fmm0_;


    // Protected Member Functions

        //- Update sub-grid eddy-viscosity
        void correctNut(const tmp<volTensorField>& gradU);

        virtual void correctNut();


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;
    typedef typename BasicMomentumTransportModel::transportModel transportModel;

    //- Runtime type information
    TypeName("dynamicLagrangian");


    // Constructors

        //- Construct from components
        dynamicLagrangian
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        dynamicLagrangian(const dynamicLagrangian&) = delete;


    //- Destructor
    virtual ~dynamicLagrangian()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        tmp<volScalarField> k(const tmp<volTensorField>& gradU) const
        {
            return
                pow(2.0*flm_/fmm_, 2.0/3.0)
              * pow(this->Ce_, -2.0/3.0)
              * sqr(this->delta())*magSqr(dev(symm(gradU)));
        }

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k(fvc::grad(this->U_));
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                this->nut_ + this->nu()
            );
        }

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const dynamicLagrangian&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicLagrangian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
