/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::SpecificIncompressibleTurbulenceModel

Description
    Templated abstract base class for specific (RAS/LES) incompressible
    turbulence models

SourceFiles
    SpecificIncompressibleTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef SpecificIncompressibleTurbulenceModel_H
#define SpecificIncompressibleTurbulenceModel_H

#include "IncompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class SpecificIncompressibleTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicIncompressibleTurbulenceModel>
class SpecificIncompressibleTurbulenceModel
:
    public BasicIncompressibleTurbulenceModel
{

public:

    typedef typename BasicIncompressibleTurbulenceModel::transportModel
        transportModel;

    // Constructors

        //- Construct from components
        SpecificIncompressibleTurbulenceModel
        (
            const word& type,
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    //- Destructor
    virtual ~SpecificIncompressibleTurbulenceModel()
    {}


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<SpecificIncompressibleTurbulenceModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SpecificIncompressibleTurbulenceModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
