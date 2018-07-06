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
    Foam::LESModels::DeardorffDiffStress

Description
    Differential SGS Stress Equation Model for incompressible and
    compressible flows

    Reference:
    \verbatim
        Deardorff, J. W. (1973).
        The use of subgrid transport equations in a three-dimensional model
        of atmospheric turbulence.
        Journal of Fluids Engineering, 95(3), 429-438.
    \endverbatim

    This SGS model uses a full balance equation for the SGS stress tensor to
    simulate the behaviour of B.

    This implementation is as described in the above paper except that the
    triple correlation model of Donaldson is replaced with the generalized
    gradient diffusion model of Daly and Harlow:
    \verbatim
        Daly, B. J., & Harlow, F. H. (1970).
        Transport equations in turbulence.
        Physics of Fluids (1958-1988), 13(11), 2634-2649.
    \endverbatim
    with the default value for the coefficient Cs of 0.25 from
    \verbatim
        Launder, B. E., Reece, G. J., & Rodi, W. (1975).
        Progress in the development of a Reynolds-stress turbulence closure.
        Journal of fluid mechanics, 68(03), 537-566.
    \endverbatim

SourceFiles
    DeardorffDiffStress.C

\*---------------------------------------------------------------------------*/

#ifndef DeardorffDiffStress_H
#define DeardorffDiffStress_H

#include "LESModel.H"
#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class DeardorffDiffStress Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class DeardorffDiffStress
:
    public ReynoldsStress<LESModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        DeardorffDiffStress(const DeardorffDiffStress&);
        void operator=(const DeardorffDiffStress&);


protected:

    // Protected data

        // Model constants

            dimensionedScalar Ck_;
            dimensionedScalar Cm_;
            dimensionedScalar Ce_;
            dimensionedScalar Cs_;


    // Protected Member Functions

        //- Update the eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("DeardorffDiffStress");


    // Constructors

        //- Constructor from components
        DeardorffDiffStress
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
    virtual ~DeardorffDiffStress()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Correct sub-grid stress, eddy-Viscosity and related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DeardorffDiffStress.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
