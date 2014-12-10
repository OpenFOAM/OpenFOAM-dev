/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::RASModels::buoyantKEpsilon

Group
    grpRASTurbulence

Description
    k-epsilon model for the gas-phase in a two-phase system
    supporting phase-inversion.

    In the limit that the gas-phase fraction approaches zero a contribution from
    the other phase is blended into the k and epsilon equations up to the
    phase-fraction of alphaInversion at which point phase-inversion is
    considered to have occurred and the model reverts to the pure single-phase
    form.

    This model is unpublished and is provided as a stable numerical framework
    on which a more physical model may be built.

    The default model coefficients correspond to the following:
    \verbatim
        buoyantKEpsilonCoeffs
        {
            Cmu             0.09;
            C1              1.44;
            C2              1.92;
            C3              -0.33;
            Cg              0.85;
            sigmak          1.0;
            sigmaEps        1.3;
            alphaInversion  0.7;
        }
    \endverbatim

SourceFiles
    buoyantKEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef buoyantKEpsilon_H
#define buoyantKEpsilon_H

#include "kEpsilon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class buoyantKEpsilon Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class buoyantKEpsilon
:
    public kEpsilon<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        buoyantKEpsilon(const buoyantKEpsilon&);
        buoyantKEpsilon& operator=(const buoyantKEpsilon&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cg_;

    // Protected Member Functions

        tmp<volScalarField> Gcoef() const;

        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("buoyantKEpsilon");


    // Constructors

        //- Construct from components
        buoyantKEpsilon
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
    virtual ~buoyantKEpsilon()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "buoyantKEpsilon.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
