/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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
    Foam::solvers::XiFluid

Description
    Solver module for compressible premixed/partially-premixed combustion with
    turbulence modelling.

    Combusting RANS code using the b-Xi two-equation model.
    Xi may be obtained by either the solution of the Xi transport
    equation or from an algebraic expression.

    Reference:
    \verbatim
        Weller, H. G. (1993).
        The development of a new flame area combustion model
        using conditional averaging.
        Thermo-fluids section report TF 9307.
    \endverbatim

    Both approaches are based on Gulder's flame speed correlation which has been
    shown to be appropriate by comparison with the results from the spectral
    model.

    Reference:
    \verbatim
        Weller, H. G., Marooney, C. J., & Gosman, A. D. (1991, January).
        A new spectral method for calculation of the time-varying area
        of a laminar flame in homogeneous turbulence.
        In Symposium (International) on Combustion
        (Vol. 23, No. 1, pp. 629-636). Elsevier.
    \endverbatim

    Strain effects are incorporated directly into the Xi equation
    but not in the algebraic approximation.  Further work need to be
    done on this issue, particularly regarding the enhanced removal rate
    caused by flame compression.  Analysis using results of the spectral
    model will be required.

    For cases involving very lean Propane flames or other flames which are
    very strain-sensitive, a transport equation for the laminar flame
    speed is present.  This equation is derived using heuristic arguments
    involving the strain time scale and the strain-rate at extinction.
    the transport velocity is the same as that for the Xi equation.

    Reference:
    \verbatim
        Weller, H. G., Tabor, G., Gosman, A. D., & Fureby, C. (1998, January).
        Application of a flame-wrinkling LES combustion model
        to a turbulent mixing layer.
        In Symposium (International) on combustion
        (Vol. 27, No. 1, pp. 899-907). Elsevier.
    \endverbatim

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Reference:
    \verbatim
        Greenshields, C. J., & Weller, H. G. (2022).
        Notes on Computational Fluid Dynamics: General Principles.
        CFD Direct Ltd.: Reading, UK.
    \endverbatim

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, chemical reactions,
    combustion, Lagrangian particles, radiation, surface film etc. and
    constraining or limiting the solution.

SourceFiles
    XiFluid.C

See also
    Foam::solvers::fluidSolver
    Foam::solvers::isothermalFluid

\*---------------------------------------------------------------------------*/

#ifndef XiFluid_H
#define XiFluid_H

#include "isothermalFluid.H"
#include "psiuMulticomponentThermo.H"
#include "RASThermophysicalTransportModel.H"
#include "unityLewisEddyDiffusivity.H"
#include "SuModel.H"
#include "XiModel.H"
#include "ignition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class XiFluid Declaration
\*---------------------------------------------------------------------------*/

class XiFluid
:
    public isothermalFluid
{

protected:

    // Thermophysical properties

        psiuMulticomponentThermo& thermo_;


    // Composition

        //- Reference to the combustion regress variable
        volScalarField& b_;

        //- Set of fields used for the multivariate convection scheme
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;


    // Thermophysical transport

        turbulenceThermophysicalTransportModels::unityLewisEddyDiffusivity
        <
            RASThermophysicalTransportModel
            <
                ThermophysicalTransportModel
                <
                    compressibleMomentumTransportModel,
                    fluidThermo
                >
            >
        > thermophysicalTransport;


    // Reactions

        //- Dictionary of combustion model coefficients
        IOdictionary combustionProperties;

        //- Strained laminar flame-speed model
        autoPtr<SuModel> SuModel_;

        //- Flame wrinkling model
        autoPtr<XiModel> XiModel_;

        //- Ignition model
        ignition ign;


    // Protected member functions

        //- Solve the ft equation for partially-premixed mixtures
        void ftSolve
        (
            const fv::convectionScheme<scalar>& mvConvection,
            const volScalarField& Db
        );

        //- Calculate and return the turbulent flame-speed kernel correction
        dimensionedScalar StCorr
        (
            const volScalarField& c,
            const surfaceScalarField& nf,
            const dimensionedScalar& dMgb
        ) const;

        //- Solve the Xi and regress variable equations
        void bSolve
        (
            const fv::convectionScheme<scalar>& mvConvection,
            const volScalarField& Db
        );

        //- Solve the unburnt energy equation
        void EauSolve
        (
            const fv::convectionScheme<scalar>& mvConvection,
            const volScalarField& Db
        );

        //- Solve the energy equation
        void EaSolve
        (
            const fv::convectionScheme<scalar>& mvConvection,
            const volScalarField& Db
        );


public:

        //- Reference to the fluid thermophysical properties
        const psiuMulticomponentThermo& thermo;

        //- Reference to the combustion regress variable
        //  obtained from the combustion mixture
        const volScalarField& b;

        //- Laminar flame speed
        const volScalarField& Su;

        //- Flame wrinkling
        const volScalarField& Xi;


    //- Runtime type information
    TypeName("XiFluid");


    // Constructors

        //- Construct from region mesh
        XiFluid(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        XiFluid(const XiFluid&) = delete;


    //- Destructor
    virtual ~XiFluid();


    // Member Functions

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const XiFluid&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
