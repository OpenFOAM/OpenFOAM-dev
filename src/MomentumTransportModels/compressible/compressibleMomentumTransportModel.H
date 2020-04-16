/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Abstract base class for turbulence models (RAS, LES and laminar).

SourceFiles
    compressibleMomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMomentumTransportModel_H
#define compressibleMomentumTransportModel_H

#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

// Transport properties traits class to provide the thermo model
template<class TransportType>
class TransportTraits
{
public:

    typedef TransportType thermoModel;

    static const thermoModel& thermo(const TransportType& tt)
    {
        return tt;
    }
};


/*---------------------------------------------------------------------------*\
             Class compressibleMomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class compressibleMomentumTransportModel
:
    public momentumTransportModel
{

protected:

    // Protected data

        const volScalarField& rho_;


public:

    //- Runtime type information
    TypeName(momentumTransportModel::typeName_());


    // Constructors

        //- Construct from components
        compressibleMomentumTransportModel
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        compressibleMomentumTransportModel
        (
            const compressibleMomentumTransportModel&
        ) = delete;


    //- Destructor
    virtual ~compressibleMomentumTransportModel()
    {}


    // Member Functions

        //- Return the density field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return the volumetric flux field
        virtual tmp<surfaceScalarField> phi() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devTau() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const compressibleMomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
