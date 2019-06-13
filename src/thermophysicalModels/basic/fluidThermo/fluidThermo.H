/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::fluidThermo

Description
    Fundamental fluid thermodynamic properties

SourceFiles
    fluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef fluidThermo_H
#define fluidThermo_H

#include "basicThermo.H"
#include "compressibleTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class fluidThermo Declaration
\*---------------------------------------------------------------------------*/

class fluidThermo
:
    public basicThermo,
    public compressibleTransportModel
{

public:

    //- Runtime type information
    TypeName("fluidThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fluidThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );

    // Constructors

        //- Construct from mesh and phase name
        fluidThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and phase name
        fluidThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );


    //- Selector
    static autoPtr<fluidThermo> New
    (
        const fvMesh&,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~fluidThermo();


    // Member Functions

        // Access to thermodynamic state variables

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual void correctRho(const volScalarField& deltaRho) = 0;

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const = 0;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
