/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::thermoCloud

Description
    Virtual abstract base class for templated ThermoCloud

SourceFiles
    thermoCloud.C

\*---------------------------------------------------------------------------*/

#ifndef thermoCloud_H
#define thermoCloud_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class thermoCloud Declaration
\*---------------------------------------------------------------------------*/

class thermoCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermoCloud(const thermoCloud&);

        //- Disallow default bitwise assignment
        void operator=(const thermoCloud&);


public:

    //- Runtime type information
    TypeName("thermoCloud");

    // Constructors

        //- Null constructor
        thermoCloud();


    //- Destructor
    virtual ~thermoCloud();


    // Member Functions

        // Access

            // Radiation

                //- Return tmp equivalent particulate emission
                virtual tmp<volScalarField> Ep() const = 0;

                //- Return tmp equivalent particulate absorption
                virtual tmp<volScalarField> ap() const = 0;

                //- Return tmp equivalent particulate scattering factor
                virtual tmp<volScalarField> sigmap() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
