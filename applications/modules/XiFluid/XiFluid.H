/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2026 OpenFOAM Foundation
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

    Combusting RANS code using the Weller b-Xi two-equation combustion model.
    Xi may be obtained by either the solution of the Xi transport equation or
    from an algebraic expression.

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

    For inhomogeneous mixtures, in addition to the regress variable \c b, it is
    necessary to solve for the mixture-fraction \c ft provided by the \c
    leanInhomogeneousMixture and also the fuel mass-fraction \c fu if there are
    rich regions in the mixture, provided by the \c inhomogeneousMixture.
    Details of the extension of the Weller b-Xi combustion model to non-premixed
    combustion can be found in the Technical Report TR/HGW/03.

    Reference:
    \verbatim
        Weller, H. G. (2002, August).
        The Application of the Weller Combustion Models to
        Non-Premixed Combustion.
        (Technical Report TR/HGW/03)
    \endverbatim

    For inhomogeneous mixtures with exhaust gas re-circulation it is necessary
    to additionally solve for the recirculated exhaust gas mass-fraction \c
    egr which is provided by the \c inhomogeneousEGRMixture mixture.

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
    Foam::solvers::basicFluidSolver
    Foam::solvers::isothermalFluid

\*---------------------------------------------------------------------------*/

#ifndef XiFluid_H
#define XiFluid_H

#include "isothermalFluid.H"
#include "ubRhoThermo.H"
#include "reactionModel.H"
#include "ubMomentumTransportModel.H"
#include "PhaseThermophysicalTransportModel.H"
#include "SuModel.H"
#include "XiModel.H"

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

        ubRhoThermo& thermo_;


    // Unburnt and burnt gas momentumTransport

        RASModels::ubMomentumTransportModel uMomentumTransport_;
        RASModels::ubMomentumTransportModel bMomentumTransport_;


    // Unburnt and burnt gas Thermophysical transport

        typedef PhaseThermophysicalTransportModel
        <
            phaseCompressible::momentumTransportModel,
            fluidThermo
        > thermophysicalTransportModel;

        autoPtr<thermophysicalTransportModel> uThermophysicalTransport_;
        autoPtr<thermophysicalTransportModel> bThermophysicalTransport_;


    // Unburnt and burnt gas reactions

        autoPtr<reactionModel> uReaction_;
        autoPtr<reactionModel> bReaction_;


    // Reactions

        //- Dictionary of combustion model coefficients
        IOdictionary combustionProperties;

        bool ignited_;

        //- Minimum b below which a numerical compensation burn-out term
        //  is applied to ensure b -> 0 and the flame speed is correct.
        //  Defaults to 1e-3
        scalar bMin_;

        //- |grad(b)| division stabilisation coefficient.
        //  Defaults to 1e-3
        scalar mgbCoeff_;

        //- Strained laminar flame-speed model
        autoPtr<SuModel> SuModel_;

        //- Flame wrinkling model
        autoPtr<XiModel> XiModel_;


    // Protected member functions

        //- Solve the regress variable and Xi equations
        //  and unburnt and burnt gas properties
        void burn();

        //- Solve for the unburnt gas thermophysical state
        void uSolve
        (
            const volScalarField::Internal& bStab,
            const surfaceScalarField& phib,
            const volScalarField::Internal& bSource
        );

        //- Solve for the unburnt gas thermophysical state
        void bSolve
        (
            const volScalarField::Internal& cStab,
            const surfaceScalarField& phic,
            const volScalarField::Internal& bSource
        );

        // Advective-diffusive stabilisation for b,c -> 0
        tmp<fvScalarMatrix> fvmStab
        (
            const volScalarField& bc,
            const volScalarField::Internal& bcStab,
            const volScalarField& D,
            volScalarField& f
        );

        //- Solve for an unburnt or burnt gas property
        void ubSolve
        (
            volScalarField& f,
            const word& fName,
            const volScalarField& alpha,
            const volScalarField& bc,
            const volScalarField::Internal& bcStab,
            const surfaceScalarField& phibc,
            const volScalarField& D,
            const thermophysicalTransportModel& thermophysicalTransport,
            const fvScalarMatrix& source,
            const reactionModel& reaction
        );

        //- Solve for a unburnt gas specie
        void uSolve
        (
            volScalarField& fu,
            const word& fuName,
            const volScalarField::Internal& bStab,
            const surfaceScalarField& phib,
            const fvScalarMatrix& source
        );

        //- Solve for a burnt gas specie
        void bSolve
        (
            volScalarField& fb,
            const word& fbName,
            const volScalarField::Internal& cStab,
            const surfaceScalarField& phic,
            const fvScalarMatrix& source
        );

        //- Solve the unburnt enthalpy equation
        void HuSolve
        (
            const volScalarField::Internal& bStab,
            const surfaceScalarField& phib,
            const volScalarField::Internal& bSource
        );

        //- Solve the burnt enthalpy equation
        void HbSolve
        (
            const volScalarField::Internal& cStab,
            const surfaceScalarField& phic,
            const volScalarField::Internal& bSource
        );


public:

        //- Reference to the fluid thermophysical properties
        const ubRhoThermo& thermo;

        const compressible::momentumTransportModel& momentumTransport;

        const thermophysicalTransportModel& uThermophysicalTransport;
        const thermophysicalTransportModel& bThermophysicalTransport;

        //- Reference to the combustion regress variable
        //  obtained from the combustion mixture
        const volScalarField& b;

        //- Reference to the burnt gas thermophysical properties
        const uRhoMulticomponentThermo& uThermo;

        //- Reference to the combustion progress variable
        const volScalarField& c;

        //- Reference to the burnt gas thermophysical properties
        const bRhoMulticomponentThermo& bThermo;

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

        //- Minimum b below which a numerical compensation burn-out term
        //  is applied to ensure b -> 0 and the flame speed is correct
        scalar bMin() const
        {
            return bMin_;
        }

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Correct the thermophysical transport
        virtual void thermophysicalTransportCorrector();

        //- Reset b-Xi and thermodynamics to the unburnt state
        //  for multi-cycle simulations
        virtual void reset();


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
