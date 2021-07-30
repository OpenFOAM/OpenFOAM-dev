/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::laminarModels::generalisedNewtonianViscosityModels::Newtonian

Description
    Newtonian viscosity model which returns the fluid Newtonian viscosity

SourceFiles
    NewtonianViscosityModel.C

\*---------------------------------------------------------------------------*/

#ifndef NewtonianViscosityModel_H
#define NewtonianViscosityModel_H

#include "generalisedNewtonianViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalisedNewtonianViscosityModels
{

/*---------------------------------------------------------------------------*\
                Class Newtonian Declaration
\*---------------------------------------------------------------------------*/

class Newtonian
:
    public generalisedNewtonianViscosityModel
{

public:

    //- Runtime type information
    TypeName("Newtonian");


    // Constructors

        //- Construct from components
        Newtonian
        (
            const dictionary& viscosityProperties,
            const viscosity& viscosity,
            const volVectorField& U
        );

        //- Disallow default bitwise copy construction
        Newtonian
        (
            const Newtonian&
        ) = delete;


    //- Destructor
    virtual ~Newtonian()
    {}


    // Member Functions

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return viscosity_.nu();
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return viscosity_.nu(patchi);
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const Newtonian&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace generalisedNewtonianViscosityModels
} // End namespace laminarModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
