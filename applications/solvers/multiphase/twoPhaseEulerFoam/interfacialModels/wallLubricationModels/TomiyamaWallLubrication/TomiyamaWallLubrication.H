/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::wallLubricationModels::TomiyamaWallLubrication

Description
    Wall lubrication model of Tomiyama.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
        in a Gas-Liquid System within an Euler-Euler Framework"
        M Otromke
        PhD Thesis
        April 2013
    \endverbatim

    \verbatim
        "Struggle with Computational Bubble Dynamics"
        A Tomiyama
        Multiphase Science and Technology
        Volume 10, Issue 4, Pages 369-405, 1998
    \endverbatim

SourceFiles
    TomiyamaWallLubrication.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaWallLubrication_H
#define TomiyamaWallLubrication_H

#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace wallLubricationModels
{

/*---------------------------------------------------------------------------*\
                            Class TomiyamaWallLubrication Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaWallLubrication
:
    public wallLubricationModel
{
private:

    // Private data

        //- Characteristic channel dimension
        const dimensionedScalar D_;


public:

    //- Runtime type information
    TypeName("Tomiyama");


    // Constructors

        //- Construct from components
        TomiyamaWallLubrication
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~TomiyamaWallLubrication();


    // Member Functions

        //- Wall lubrication force
        tmp<volVectorField> F() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
