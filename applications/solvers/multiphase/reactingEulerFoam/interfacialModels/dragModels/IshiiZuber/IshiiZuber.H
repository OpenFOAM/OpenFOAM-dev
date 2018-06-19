/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::dragModels::IshiiZuber

Description
    Ishii and Zuber (1979) drag model for dense dispersed bubbly flows.

    Reference:
    \verbatim
        "Drag Coefficient and relative velocity in bubbly, droplet and
         particulate flows",
        Ishii, M., Zuber, N.,
        AIChE Journal 5, Vol. 25, 1979, pp. 843-855.
    \endverbatim

SourceFiles
    IshiiZuber.C

\*---------------------------------------------------------------------------*/

#ifndef IshiiZuber_H
#define IshiiZuber_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class IshiiZuber Declaration
\*---------------------------------------------------------------------------*/

class IshiiZuber
:
    public dragModel
{

public:

    //- Runtime type information
    TypeName("IshiiZuber");


    // Constructors

        //- Construct from a dictionary and a phase pair
        IshiiZuber
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~IshiiZuber();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
