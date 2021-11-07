/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::compressible::RASModels::PDRkEpsilon

Description
    Standard k-epsilon turbulence model with additional source terms
    corresponding to PDR basic drag model (\link basic.H \endlink)

    The default model coefficients correspond to the following:
    @verbatim
        PDRkEpsilonCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            C3          -0.33;  // only for compressible
            C4          0.1;
            sigmak      1.0;    // only for compressible
            sigmaEps    1.3;
            Prt         1.0;    // only for compressible
        }
    @endverbatim

    The turbulence source term \f$ G_{R} \f$ appears in the
    \f$ \kappa-\epsilon \f$ equation for the generation of turbulence due to
    interaction with unresolved obstacles.

    In the \f$ \epsilon  \f$ equation \f$ C_{1} G_{R} \f$ is added as a source
    term.

    In the \f$ \kappa \f$ equation \f$ G_{R} \f$ is added as a source term.

SourceFiles
    PDRkEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef compressible_PDRkEpsilon_H
#define compressible_PDRkEpsilon_H

#include "kEpsilon.H"
#include "compressibleMomentumTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                         Class PDRkEpsilon Declaration
\*---------------------------------------------------------------------------*/

class PDRkEpsilon
:
    public Foam::RASModels::kEpsilon<compressible::momentumTransportModel>
{
    // Private Data

        // Model coefficients
        dimensionedScalar C4_;


public:

    //- Runtime type information
    TypeName("PDRkEpsilon");


    // Constructors

        //- Construct from components
        PDRkEpsilon
        (
            const geometricOneField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& transport,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~PDRkEpsilon();


    // Member Functions

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Read momentumTransport dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
