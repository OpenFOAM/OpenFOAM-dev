/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::eddyViscosity

Description
    Eddy viscosity turbulence model base class

SourceFiles
    eddyViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef eddyViscosity_H
#define eddyViscosity_H

#include "linearViscousStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class eddyViscosity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class eddyViscosity
:
    public linearViscousStress<BasicTurbulenceModel>
{

protected:

    // Protected data

        // Fields

            volScalarField nut_;


    // Protected Member Functions

        virtual void correctNut() = 0;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct from components
        eddyViscosity
        (
            const word& modelName,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    //- Destructor
    virtual ~eddyViscosity()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read() = 0;

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the turbulence viscosity on patch
        virtual tmp<scalarField> nut(const label patchi) const
        {
            return nut_.boundaryField()[patchi];
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const = 0;

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Validate the turbulence fields after construction
        //  Update turbulence viscosity and other derived fields as requires
        virtual void validate();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "eddyViscosity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
