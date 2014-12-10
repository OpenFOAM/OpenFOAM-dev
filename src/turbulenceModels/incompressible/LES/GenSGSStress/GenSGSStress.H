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
    Foam::incompressible::LESModels::GenSGSStress

Group
    grpIcoLESTurbulence

Description
    General base class for all incompressible models that directly
    solve for the SGS stress tensor B.

    Contains tensor fields B (the SGS stress tensor) as well as scalar
    fields for k (SGS turbulent energy) gamma (SGS viscosity) and epsilon
    (SGS dissipation).

SourceFiles
    GenSGSStress.C

\*---------------------------------------------------------------------------*/

#ifndef GenSGSStress_H
#define GenSGSStress_H

#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class GenSGSStress Declaration
\*---------------------------------------------------------------------------*/

class GenSGSStress
:
    virtual public LESModel
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        GenSGSStress(const GenSGSStress&);
        GenSGSStress& operator=(const GenSGSStress&);


protected:

        dimensionedScalar ce_;

        dimensionedScalar couplingFactor_;

        volSymmTensorField B_;
        volScalarField nuSgs_;


public:

    //- Partial Runtime type information
    static const word typeName;

    // Constructors

        //- Construct from components
        GenSGSStress
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~GenSGSStress()
    {}


    // Member Functions

        //- Return the SGS turbulent kinetic energy.
        virtual tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (   
                        "k",
                        runTime_.timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    0.5*tr(B_)
                )
            );
        }

        //- Return the SGS turbulent dissipation.
        virtual tmp<volScalarField> epsilon() const
        {
            const volScalarField K(k());
            
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
                    ce_*K*sqrt(K)/delta()
                )
            );
        }

        //- Return the SGS viscosity.
        virtual tmp<volScalarField> nuSgs() const
        {
            return nuSgs_;
        }

        //- Return the sub-grid stress tensor.
        virtual tmp<volSymmTensorField> B() const
        {
            return B_;
        }

        //- Return the effective sub-grid turbulence stress tensor
        //  including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the deviatoric part of the effective sub-grid
        //  turbulence stress tensor including the laminar stress
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Read LESProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
