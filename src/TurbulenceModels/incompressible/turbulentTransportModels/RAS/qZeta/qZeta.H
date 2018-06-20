/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::incompressible::RASModels::qZeta

Description
    Gibson and Dafa'Alla's q-zeta two-equation low-Re turbulence model
    for incompressible flows

    This turbulence model is described in:
    \verbatim
        Dafa'Alla, A.A., Juntasaro, E. & Gibson, M.M. (1996).
        Calculation of oscillating boundary layers with the
        q-zeta turbulence model.
        Engineering Turbulence Modelling and Experiments 3:
        Proceedings of the Third International Symposium,
        Crete, Greece, May 27-29, 141.
    \endverbatim
    which is a development of the original q-zeta model described in:
    \verbatim
        Gibson, M. M., & Dafa'Alla, A. A. (1995).
        Two-equation model for turbulent wall flow.
        AIAA journal, 33(8), 1514-1518.
    \endverbatim

SourceFiles
    qZeta.C

\*---------------------------------------------------------------------------*/

#ifndef qZeta_H
#define qZeta_H

#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class qZeta Declaration
\*---------------------------------------------------------------------------*/

class qZeta
:
    public eddyViscosity<incompressible::RASModel>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaZeta_;
            Switch anisotropic_;

        //- Lower limit of q
        dimensionedScalar qMin_;

        //- Lower limit of zeta
        dimensionedScalar zetaMin_;

        // Fields

            volScalarField k_;
            volScalarField epsilon_;

            volScalarField q_;
            volScalarField zeta_;


    // Protected Member Functions

        tmp<volScalarField> fMu() const;
        tmp<volScalarField> f2() const;
        virtual void correctNut();


public:

    //- Runtime type information
    TypeName("qZeta");

    // Constructors

        //- Construct from components
        qZeta
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~qZeta()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the lower allowable limit for q (default: small)
        const dimensionedScalar& qMin() const
        {
            return qMin_;
        }

        //- Return the lower allowable limit for zeta (default: small)
        const dimensionedScalar& zetaMin() const
        {
            return zetaMin_;
        }

        //- Allow qMin to be changed
        dimensionedScalar& qMin()
        {
            return qMin_;
        }

        //- Allow zetaMin to be changed
        dimensionedScalar& zetaMin()
        {
            return zetaMin_;
        }

        //- Return the effective diffusivity for q
        tmp<volScalarField> DqEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DqEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DzetaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DzetaEff", nut_/sigmaZeta_ + nu())
            );
        }

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

        virtual const volScalarField& q() const
        {
            return q_;
        }

        virtual const volScalarField& zeta() const
        {
            return zeta_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
