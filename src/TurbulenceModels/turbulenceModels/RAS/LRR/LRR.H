/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::RASModels::LRR

Description
    Launder, Reece and Rodi Reynolds-stress turbulence model for
    incompressible and compressible flows.

    Reference:
    \verbatim
        Launder, B. E., Reece, G. J., & Rodi, W. (1975).
        Progress in the development of a Reynolds-stress turbulence closure.
        Journal of fluid mechanics, 68(03), 537-566.
    \endverbatim

    Including the recommended generalized gradient diffusion model of
    Daly and Harlow:
    \verbatim
        Daly, B. J., & Harlow, F. H. (1970).
        Transport equations in turbulence.
        Physics of Fluids (1958-1988), 13(11), 2634-2649.
    \endverbatim

    Optional Gibson-Launder wall-reflection is also provided:
    \verbatim
        Gibson, M. M., & Launder, B. E. (1978).
        Ground effects on pressure fluctuations in the
        atmospheric boundary layer.
        Journal of Fluid Mechanics, 86(03), 491-511.
    \endverbatim

    The default model coefficients are:
    \verbatim
        LRRCoeffs
        {
            Cmu             0.09;
            C1              1.8;
            C2              0.6;
            Ceps1           1.44;
            Ceps2           1.92;
            Cs              0.25;
            Ceps            0.15;

            wallReflection  yes;
            kappa           0.41
            Cref1           0.5;
            Cref2           0.3;

            couplingFactor  0.0;
        }
    \endverbatim

SourceFiles
    LRR.C

\*---------------------------------------------------------------------------*/

#ifndef LRR_H
#define LRR_H

#include "RASModel.H"
#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class LRR Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class LRR
:
    public ReynoldsStress<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        LRR(const LRR&);
        void operator=(const LRR&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;

            dimensionedScalar C1_;
            dimensionedScalar C2_;

            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar Cs_;
            dimensionedScalar Ceps_;


        // Wall-refection coefficients

            Switch wallReflection_;
            dimensionedScalar kappa_;
            dimensionedScalar Cref1_;
            dimensionedScalar Cref2_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;


    // Protected Member Functions

        //- Update the eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("LRR");


    // Constructors

        //- Construct from components
        LRR
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
    virtual ~LRR()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the effective diffusivity for R
        tmp<volSymmTensorField> DREff() const;

        //- Return the effective diffusivity for epsilon
        tmp<volSymmTensorField> DepsilonEff() const;

        //- Solve the turbulence equations and correct eddy-Viscosity and
        //  related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LRR.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
