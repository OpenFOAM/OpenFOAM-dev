/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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
    Foam::solidPressureThermo

Description
    Fundamental solid thermodynamic properties including pressure.

SourceFiles
    solidPressureThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidPressureThermo_H
#define solidPressureThermo_H

#include "solidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solidPressureThermo Declaration
\*---------------------------------------------------------------------------*/

class solidPressureThermo
:
    public solidThermo
{

protected:

        // Fields

            //- Pressure [Pa]
            volScalarField& p_;


public:

    //- Runtime type information
    TypeName("solidPressureThermo");

    // Constructors

        //- Construct from mesh and phase name
        solidPressureThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh, dictionary and phase name
        solidPressureThermo
        (
            const fvMesh&,
            const dictionary& dict,
            const word& phaseName
        );


    //- Destructor
    virtual ~solidPressureThermo();


    // Member Functions

        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            virtual volScalarField& p();

            //- Pressure [Pa]
            virtual const volScalarField& p() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
