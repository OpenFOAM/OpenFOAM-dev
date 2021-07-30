/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
    Foam::laminarModels::generalisedNewtonianViscosityModels

Description
    A namespace for the generalised Newtonian viscosity model implementations.

Class
    Foam::laminarModels::generalisedNewtonianViscosityModel

Description
    An abstract base class for generalised Newtonian viscosity models.

SourceFiles
    generalisedNewtonianViscosityModel.C
    generalisedNewtonianViscosityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef generalisedNewtonianViscosityModel_H
#define generalisedNewtonianViscosityModel_H

#include "viscosity.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{

/*---------------------------------------------------------------------------*\
             Class generalisedNewtonianViscosityModel Declaration
\*---------------------------------------------------------------------------*/

class generalisedNewtonianViscosityModel
:
    public viscosity
{
protected:

        //- Reference to the fluid Newtonian viscosity
        const viscosity& viscosity_;

        //- Reference to the velocity field
        const volVectorField& U_;


public:

    //- Runtime type information
    TypeName("generalisedNewtonianViscosityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            generalisedNewtonianViscosityModel,
            dictionary,
            (
                const dictionary& viscosityProperties,
                const viscosity& viscosity,
                const volVectorField& U
            ),
            (viscosityProperties, viscosity, U)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<generalisedNewtonianViscosityModel> New
        (
            const dictionary& viscosityProperties,
            const viscosity& viscosity,
            const volVectorField& U
        );


    // Constructors

        //- Construct from components
        generalisedNewtonianViscosityModel
        (
            const dictionary& viscosityProperties,
            const viscosity& viscosity,
            const volVectorField& U
        );

        //- Disallow default bitwise copy construction
        generalisedNewtonianViscosityModel
        (
            const generalisedNewtonianViscosityModel&
        ) = delete;


    //- Destructor
    virtual ~generalisedNewtonianViscosityModel()
    {}


    // Member Functions

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties) = 0;

        //- Correct the laminar viscosity
        virtual void correct() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const generalisedNewtonianViscosityModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
