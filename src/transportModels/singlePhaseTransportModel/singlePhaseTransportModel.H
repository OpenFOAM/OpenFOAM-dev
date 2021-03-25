/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::singlePhaseTransportModel

Description
    A simple single-phase transport model based on viscosityModel.

    Used by the incompressible single-phase solvers like simpleFoam,
    turbFoam etc.

SourceFiles
    singlePhaseTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseTransportModel_H
#define singlePhaseTransportModel_H

#include "kinematicTransportModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class viscosityModel;

/*---------------------------------------------------------------------------*\
                  Class singlePhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseTransportModel
:
    public IOdictionary,
    public kinematicTransportModel
{
    // Private Data

        autoPtr<viscosityModel> viscosityModelPtr_;


public:

    //- Runtime type information
    TypeName("singlePhaseTransportModel");


    // Constructors

        //- Construct from components
        singlePhaseTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        singlePhaseTransportModel(const singlePhaseTransportModel&);


    //- Destructor
    virtual ~singlePhaseTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const singlePhaseTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
