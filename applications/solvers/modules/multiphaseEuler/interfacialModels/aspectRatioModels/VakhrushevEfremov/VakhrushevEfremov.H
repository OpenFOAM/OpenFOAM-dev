/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
    Foam::aspectRatioModels::VakhrushevEfremov

Description
    Aspect ratio model of Vakhrushev and Efremov.

    Reference:
    \verbatim
        Vakhrushev, I. A., & Efremov, G. I. (1970).
        Interpolation formula for computing the velocities of single gas bubbles
        in liquids.
        Chemistry and Technology of Fuels and Oils, 6(5), 376-379.
    \endverbatim

SourceFiles
    VakhrushevEfremov.C

\*---------------------------------------------------------------------------*/

#ifndef VakhrushevEfremov_H
#define VakhrushevEfremov_H

#include "aspectRatioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aspectRatioModels
{

/*---------------------------------------------------------------------------*\
                      Class VakhrushevEfremov Declaration
\*---------------------------------------------------------------------------*/

class VakhrushevEfremov
:
    public aspectRatioModel
{
public:

    //- Runtime type information
    TypeName("VakhrushevEfremov");


    // Constructors

        //- Construct from a dictionary and an interface
        VakhrushevEfremov
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~VakhrushevEfremov();


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aspectRatioModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
