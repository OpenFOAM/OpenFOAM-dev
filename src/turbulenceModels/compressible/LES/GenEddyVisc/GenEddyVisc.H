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
    Foam::compressible::LESModels::GenEddyVisc

Group
    grpCmpLESTurbulence

Description
    General base class for all compressible models that can be implemented as
    an eddy viscosity, i.e. algebraic and one-equation models.

    Contains fields for k (SGS turbulent kinetic energy), gamma
    (modelled viscosity) and epsilon (SGS dissipation).

SourceFiles
    GenEddyVisc.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleGenEddyVisc_H
#define compressibleGenEddyVisc_H

#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                         Class GenEddyVisc Declaration
\*---------------------------------------------------------------------------*/

class GenEddyVisc
:
    virtual public LESModel
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        GenEddyVisc(const GenEddyVisc&);
        GenEddyVisc& operator=(const GenEddyVisc&);


protected:

    // Model coefficients

        dimensionedScalar ce_;
        dimensionedScalar Prt_;


    // Fields

        volScalarField muSgs_;
        volScalarField alphaSgs_;


public:

    //- Partial Runtime type information
    static const word typeName;

    // Constructors

        //- Construct from components
        GenEddyVisc
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const fluidThermo& thermoPhysicalModel,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~GenEddyVisc()
    {}


    // Member Functions

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (   
                        "epsilon",
                        runTime_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    ce_*k()*sqrt(k())/delta()
                )
            );
        }

        //- Return viscosity
        virtual tmp<volScalarField> muSgs() const
        {
            return muSgs_;
        }

        //- Return thermal diffusivity
        virtual tmp<volScalarField> alphaSgs() const
        {
            return alphaSgs_;
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoBeff() const;

        //- Returns div(rho*dev(B)).
        // This is the additional term due to the filtering of the NSE.
        virtual tmp<fvVectorMatrix> divDevRhoBeff(volVectorField& U) const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct(const tmp<volTensorField>& gradU);

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
