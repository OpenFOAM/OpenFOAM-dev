/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    Foam::compressibleMomentumTransportModel

Description
    Base class for single-phase compressible turbulence models.

SourceFiles
    compressibleMomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMomentumTransportModel_H
#define compressibleMomentumTransportModel_H

#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class compressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class compressibleMomentumTransportModel
:
    public momentumTransportModel
{

public:

    typedef geometricOneField alphaField;
    typedef volScalarField rhoField;


protected:

    // Protected data

        alphaField alpha_;
        const rhoField& rho_;


public:

    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            compressibleMomentumTransportModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const viscosity& viscosity
            ),
            (alpha, rho, U, alphaRhoPhi, phi, viscosity)
        );


    // Constructors

        //- Construct
        compressibleMomentumTransportModel
        (
            const word& type,
            const geometricOneField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<compressibleMomentumTransportModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const viscosity& viscosity
        );


    //- Destructor
    virtual ~compressibleMomentumTransportModel()
    {}


    // Member Functions

        //- Return the phase fraction field
        const alphaField& alpha() const
        {
            return alpha_;
        }

        //- Return the density field
        const rhoField& rho() const
        {
            return rho_;
        }

        //- Return the volumetric flux field
        virtual tmp<surfaceScalarField> phi() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devTau() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const = 0;
};


namespace compressible
{
    typedef compressibleMomentumTransportModel momentumTransportModel;

    template<class CompressibleMomentumTransportModel>
    autoPtr<CompressibleMomentumTransportModel> New
    (
        const volScalarField& rho,
        const volVectorField& U,
        const surfaceScalarField& phi,
        const viscosity& viscosity
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "compressibleMomentumTransportModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
