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
    A namespace for various generalised Newtonian viscosity model
    implementations.

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

#include "dictionary.H"
#include "volFieldsFwd.H"
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
{
protected:

    // Protected data

        dictionary viscosityProperties_;


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
                const dictionary& viscosityProperties
            ),
            (viscosityProperties)
        );


    // Selectors

        //- Return a reference to the selected viscosity model
        static autoPtr<generalisedNewtonianViscosityModel> New
        (
            const dictionary& viscosityProperties
        );


    // Constructors

        //- Construct from components
        generalisedNewtonianViscosityModel
        (
            const dictionary& viscosityProperties
        );

        //- Disallow default bitwise copy construction
        generalisedNewtonianViscosityModel
        (
            const generalisedNewtonianViscosityModel&
        );


    //- Destructor
    virtual ~generalisedNewtonianViscosityModel()
    {}


    // Member Functions

        //- Return the phase transport properties dictionary
        const dictionary& viscosityProperties() const
        {
            return viscosityProperties_;
        }

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu
        (
            const volScalarField& nu0,
            const volScalarField& strainRate
        ) const = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties) = 0;


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
