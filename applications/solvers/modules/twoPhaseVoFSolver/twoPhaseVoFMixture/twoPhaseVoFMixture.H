/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::twoPhaseVoFMixture

Description
    Class to represent a VoF mixture

SourceFiles
    twoPhaseVoFMixture.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseVoFMixture_H
#define twoPhaseVoFMixture_H

#include "VoFMixture.H"
#include "twoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class twoPhaseVoFMixture Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseVoFMixture
:
    public VoFMixture,
    public twoPhaseMixture
{

public:

    TypeName("twoPhaseVoFMixture");


    // Constructors

        //- Construct from components
        twoPhaseVoFMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~twoPhaseVoFMixture()
    {}


    // Member Functions

        //- Return the mixture density
        using twoPhaseMixture::rho;

        //- Correct the mixture properties
        using twoPhaseMixture::correct;

        //- Read base phaseProperties dictionary
        using twoPhaseMixture::read;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
