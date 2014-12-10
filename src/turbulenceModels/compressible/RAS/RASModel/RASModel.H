/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Namespace
    Foam::compressible::RASModels

Description
    Namespace for compressible RAS turbulence models.

Class
    Foam::compressible::RASModel

Description
    Abstract base class for turbulence models for compressible and combusting
    flows.

SourceFiles
    RASModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleRASModel_H
#define compressibleRASModel_H

#include "compressible/turbulenceModel/turbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fluidThermo.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
                           Class RASModel Declaration
\*---------------------------------------------------------------------------*/

class RASModel
:
    public turbulenceModel,
    public IOdictionary
{

protected:

    // Protected data

        //- Turbulence on/off flag
        Switch turbulence_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- Lower limit of k
        dimensionedScalar kMin_;

        //- Lower limit of epsilon
        dimensionedScalar epsilonMin_;

        //- Lower limit for omega
        dimensionedScalar omegaMin_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        RASModel(const RASModel&);

        //- Disallow default bitwise assignment
        void operator=(const RASModel&);


public:

    //- Runtime type information
    TypeName("RASModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RASModel,
            dictionary,
            (
                const volScalarField& rho,
                const volVectorField& U,
                const surfaceScalarField& phi,
                const fluidThermo& thermoPhysicalModel,
                const word& turbulenceModelName
            ),
            (rho, U, phi, thermoPhysicalModel, turbulenceModelName)
        );


    // Constructors

        //- Construct from components
        RASModel
        (
            const word& type,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<RASModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    //- Destructor
    virtual ~RASModel()
    {}


    // Member Functions

        // Access

            //- Return the lower allowable limit for k (default: SMALL)
            const dimensionedScalar& kMin() const
            {
                return kMin_;
            }

            //- Return the lower allowable limit for epsilon (default: SMALL)
            const dimensionedScalar& epsilonMin() const
            {
                return epsilonMin_;
            }

            //- Return the lower allowable limit for omega (default: SMALL)
            const dimensionedScalar& omegaMin() const
            {
                return omegaMin_;
            }

            //- Allow kMin to be changed
            dimensionedScalar& kMin()
            {
                return kMin_;
            }

            //- Allow epsilonMin to be changed
            dimensionedScalar& epsilonMin()
            {
                return epsilonMin_;
            }

            //- Allow omegaMin to be changed
            dimensionedScalar& omegaMin()
            {
                return omegaMin_;
            }

            //- Const access to the coefficients dictionary
            virtual const dictionary& coeffDict() const
            {
                return coeffDict_;
            }


        //- Return the effective viscosity
        virtual tmp<volScalarField> muEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("muEff", mut() + mu())
            );
        }

        //- Return the effective turbulent temperature diffusivity
        virtual tmp<volScalarField> kappaEff() const
        {
            return thermo().kappaEff(alphat());
        }

        //- Return the effective turbulent temperature diffusivity for a patch
        virtual tmp<scalarField> kappaEff(const label patchI) const
        {
            return thermo().kappaEff(alphat(patchI), patchI);
        }

        //- Return the effective turbulent thermal diffusivity
        virtual tmp<volScalarField> alphaEff() const
        {
            return thermo().alphaEff(alphat());
        }

        //- Return the effective turbulent thermal diffusivity for a patch
        virtual tmp<scalarField> alphaEff(const label patchI) const
        {
            return thermo().alphaEff(alphat(patchI), patchI);
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
