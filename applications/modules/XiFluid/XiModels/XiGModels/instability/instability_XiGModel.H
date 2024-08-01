/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::XiGModels::instability

Description
    Flame instability flame wrinkling generation rate model

SourceFiles
    instability.C

\*---------------------------------------------------------------------------*/

#ifndef instability_XiGModel_H
#define instability_XiGModel_H

#include "XiGModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{

/*---------------------------------------------------------------------------*\
                        Class instability Declaration
\*---------------------------------------------------------------------------*/

class instability
:
    public XiGModel
{
    // Private Data

        //- Flame instability wrinkling generation rate
        dimensionedScalar Gin_;

        //- Xi generation rate model due to all other processes
        autoPtr<XiGModel> XiGModel_;


protected:

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("instability");


    // Constructors

        //- Construct from components
        instability
        (
            const dictionary& dict,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        instability(const instability&) = delete;


    //- Destructor
    virtual ~instability();


    // Member Functions

        //- Return the flame-wrinkling generation rate
        virtual tmp<volScalarField> G() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const instability&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace XiGModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
