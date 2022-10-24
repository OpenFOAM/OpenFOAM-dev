/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::solidThermophysicalTransportModel

Description
    Abstract base class for solid thermophysical transport models

SourceFiles
    solidThermophysicalTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidThermophysicalTransportModel_H
#define solidThermophysicalTransportModel_H

#include "thermophysicalTransportModel.H"
#include "solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class solidThermophysicalTransportModel Declaration
\*---------------------------------------------------------------------------*/

class solidThermophysicalTransportModel
:
    public thermophysicalTransportModel
{
protected:

    // Protected data

        //- Reference to the solid thermophysical properties
        const solidThermo& thermo_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs(const word& type);


public:

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidThermophysicalTransportModel,
            dictionary,
            (
                const solidThermo& thermo
            ),
            (thermo)
        );


    // Constructors

        //- Construct from solid thermophysical properties
        solidThermophysicalTransportModel
        (
            const word& type,
            const solidThermo& thermo
        );


    // Selectors

        //- Return a reference to the selected thermophysical transport model
        static autoPtr<solidThermophysicalTransportModel> New
        (
            const solidThermo& thermo
        );


    //- Destructor
    virtual ~solidThermophysicalTransportModel()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read() = 0;

        //- Access function to solid thermophysical properties
        virtual const solidThermo& thermo() const
        {
            return thermo_;
        }

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Thermal conductivity [W/m/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal conductivity for patch [W/m/K]
        virtual tmp<scalarField> kappa(const label patchi) const;

        //- Effective thermal conductivity
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            return kappa();
        }

        //- Effective thermal conductivity
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return kappa(patchi);
        }

        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const = 0;

        //- Return the patch heat flux correction [W/m^2]
        //  For isotropic or patch-aligned thermal conductivity qCorr is null
        virtual tmp<scalarField> qCorr(const label patchi) const = 0;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const = 0;

        //- Solve the thermophysical transport model equations
        //  and correct the transport coefficients
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
