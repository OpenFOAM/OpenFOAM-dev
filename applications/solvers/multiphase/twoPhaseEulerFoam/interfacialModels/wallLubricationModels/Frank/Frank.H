/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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
    Foam::wallLubricationModels::Frank

Description
    Wall lubrication model of Frank.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
         in a Gas-Liquid System within an Euler-Euler Framework"
        Otromke, M.
        PhD Thesis
        April 2013
    \endverbatim

    \verbatim
        "Advances in Computational Fluid Dynamics (CFD) of 3-dimensional Gas-
         Liquid Multiphase Flows"
        Frank, T.
        NAFEMS Seminar: Simulation of Complex Flows (CFD), April 2005, pp. 1-18
    \endverbatim

SourceFiles
    Frank.C

\*---------------------------------------------------------------------------*/

#ifndef Frank_H
#define Frank_H

#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                            Class Frank Declaration
\*---------------------------------------------------------------------------*/

class Frank
:
    public wallLubricationModel
{
    // Private Data

        //- Coefficient d
        const dimensionedScalar Cwd_;

        //- Coefficient c
        const dimensionedScalar Cwc_;

        //- Power p
        const scalar p_;


public:

    //- Runtime type information
    TypeName("Frank");


    // Constructors

        //- Construct from components
        Frank
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Frank();


    // Member Functions

        //- Return phase-intensive wall lubrication force
        tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
