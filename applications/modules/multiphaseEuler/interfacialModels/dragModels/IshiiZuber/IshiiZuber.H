/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
        Ishii, M., & Zuber, N. (1979).
        Drag coefficient and relative velocity in bubbly, droplet or particulate
        flows.
        AIChE Journal, 25(5), 843-855.
    \endverbatim

SourceFiles
    IshiiZuber.C

\*---------------------------------------------------------------------------*/

#ifndef IshiiZuber_H
#define IshiiZuber_H

#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                         Class IshiiZuber Declaration
\*---------------------------------------------------------------------------*/

class IshiiZuber
:
    public dispersedDragModel
{

public:

    //- Runtime type information
    TypeName("IshiiZuber");


    // Constructors

        //- Construct from a dictionary and an interface
        IshiiZuber
        (
            const dictionary& dict,
            const phaseInterface& interface,
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
