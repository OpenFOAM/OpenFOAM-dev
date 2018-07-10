/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::LESModels::kOmegaSSTDES

Description
    Implementation of the k-omega-SST-DES turbulence model for
    incompressible and compressible flows.

    DES model described in:
    \verbatim
        Menter, F. R., Kuntz, M., and Langtry, R. (2003).
        Ten Years of Industrial Experience with the SST Turbulence Model.
        Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano,
        & M. Tummers, Begell House, Inc., 625 - 632.
    \endverbatim

    Optional support for zonal filtering based on F1 or F2 is provided as
    described in the paper.

    For further details of the implementation of the base k-omega-SST model
    see Foam::kOmegaSST.

See also
    Foam::kOmegaSST

SourceFiles
    kOmegaSST.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTDES_H
#define kOmegaSSTDES_H

#include "kOmegaSSTBase.H"
#include "LESModel.H"
#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSST Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTDES
:
    public Foam::kOmegaSST
    <
        LESeddyViscosity<BasicTurbulenceModel>,
        BasicTurbulenceModel
    >
{

protected:

    // Protected data

        // Model constants

            //- DES coefficient
            dimensionedScalar CDES_;

            //- Zonal filter choice
            //
            //  - 0: no filtering
            //  - 1: (1 - F1)
            //  - 2: (1 - F2)
            direction FSST_;


    // Protected Member Functions

        //- Return the turbulent length-scale
        tmp<volScalarField::Internal> Lt() const;

        //- The DES dissipation-rate multiplier with options zonal filtering
        //  based on either F1 or F2
        virtual tmp<volScalarField::Internal> FDES
        (
            const volScalarField::Internal& F1,
            const volScalarField::Internal& F2
        ) const;

        //- Return epsilon/k which for standard RAS is betaStar*omega
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField::Internal& F1,
            const volScalarField::Internal& F2
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTDES");


    // Constructors

        //- Construct from components
        kOmegaSSTDES
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
    virtual ~kOmegaSSTDES()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaSSTDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
