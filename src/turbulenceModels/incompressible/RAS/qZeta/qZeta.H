/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Group
    grpIcoRASTurbulence

Description
    Gibson and Dafa'Alla's q-zeta two-equation low-Re turbulence model
    for incompressible flows

    References:
    \verbatim
        "Calculation of oscillating boundary layers with the
        q-[zeta] turbulence model"
        A.A. Dafa'Alla, E. Juntasaro and M.M. Gibson

        Engineering Turbulence Modelling and Experiments 3:
        Proceedings of the Third International Symposium,
        Crete, Greece, May 27-29, 1996, p. 141.

        Editors: Wolfgang Rodi and G. Bergeles
    \endverbatim

SourceFiles
    qZeta.C

\*---------------------------------------------------------------------------*/

#ifndef qZeta_H
#define qZeta_H

#include "RASModel.H"

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
    public RASModel
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

            volScalarField nut_;


    // Protected Member Functions

        tmp<volScalarField> fMu() const;
        tmp<volScalarField> f2() const;


public:

    //- Runtime type information
    TypeName("qZeta");

    // Constructors

        //- Construct from components
        qZeta
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~qZeta()
    {}


    // Member Functions

       // Access

            //- Return the lower allowable limit for q (default: SMALL)
            const dimensionedScalar& qMin() const
            {
                return qMin_;
            }

            //- Return the lower allowable limit for zeta (default: SMALL)
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


        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
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

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        virtual const volScalarField& q() const
        {
            return q_;
        }

        virtual const volScalarField& zeta() const
        {
            return zeta_;
        }

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
