/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::SuModels::uniformConstant

Description
    Simple uniform constant laminar flame speed \c Su model

SourceFiles
    uniformConstant_SuModel.C

\*---------------------------------------------------------------------------*/

#ifndef uniformConstant_SuModel_H
#define uniformConstant_SuModel_H

#include "SuModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace SuModels
{

/*---------------------------------------------------------------------------*\
                            Class uniformConstant Declaration
\*---------------------------------------------------------------------------*/

class uniformConstant
:
    public SuModel
{
    // Private Data

        //- Specified laminar flame speed
        dimensionedScalar Su_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("uniformConstant");


    // Constructors

        //- Construct from components
        uniformConstant
        (
            const dictionary& dict,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence
        );

        //- Disallow default bitwise copy construction
        uniformConstant(const uniformConstant&) = delete;


    //- Destructor
    virtual ~uniformConstant();


    // Member Functions

        //- Correct the laminar flame speed
        virtual void correct()
        {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uniformConstant&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SuModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
