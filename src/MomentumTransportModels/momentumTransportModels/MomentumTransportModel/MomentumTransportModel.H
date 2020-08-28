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
    Foam::MomentumTransportModel

Description
    Templated abstract base class for turbulence models

SourceFiles
    MomentumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef MomentumTransportModel_H
#define MomentumTransportModel_H

#include "momentumTransportModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class MomentumTransportModel Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Alpha,
    class Rho,
    class BasicMomentumTransportModel,
    class TransportModel
>
class MomentumTransportModel
:
    public BasicMomentumTransportModel
{

public:

    typedef Alpha alphaField;
    typedef Rho rhoField;
    typedef TransportModel transportModel;


protected:

    // Protected data

        const alphaField& alpha_;
        const transportModel& transport_;


public:

    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            MomentumTransportModel,
            dictionary,
            (
                const alphaField& alpha,
                const rhoField& rho,
                const volVectorField& U,
                const surfaceScalarField& alphaRhoPhi,
                const surfaceScalarField& phi,
                const transportModel& transport
            ),
            (alpha, rho, U, alphaRhoPhi, phi, transport)
        );


    // Constructors

        //- Construct
        MomentumTransportModel
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );

        //- Disallow default bitwise copy construction
        MomentumTransportModel(const MomentumTransportModel&) = delete;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<MomentumTransportModel> New
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport
        );


    //- Destructor
    virtual ~MomentumTransportModel()
    {}


    // Member Functions

        //- Access function to phase fraction
        const alphaField& alpha() const
        {
            return alpha_;
        }

        //- Access function to incompressible transport model
        const transportModel& transport() const
        {
            return transport_;
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MomentumTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MomentumTransportModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
