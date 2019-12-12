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
    Foam::regionModels::surfaceFilmModels::primaryRadiation

Description
    Radiation model whereby the radiative heat flux is mapped from the primary
    region

SourceFiles
    primaryRadiation.C

\*---------------------------------------------------------------------------*/

#ifndef primaryRadiation_H
#define primaryRadiation_H

#include "filmRadiationModel.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                     Class primaryRadiation Declaration
\*---------------------------------------------------------------------------*/

class primaryRadiation
:
    public filmRadiationModel
{
    // Private Data

        //- Incident radiative flux mapped from the primary region [kg/s^3]
        volScalarField::Internal qinFilm_;


public:

    //- Runtime type information
    TypeName("primaryRadiation");


    // Constructors

        //- Construct from surface film model and dictionary
        primaryRadiation
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        primaryRadiation(const primaryRadiation&) = delete;


    //- Destructor
    virtual ~primaryRadiation();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct();

            //- Return the radiation sensible enthalpy source
            //  Also updates qrNet
            virtual tmp<volScalarField::Internal> Shs();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const primaryRadiation&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
