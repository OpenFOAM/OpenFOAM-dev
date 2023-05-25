/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::dragModels::Ergun

Description

    Reference:
    \verbatim
        Enwald, H., Peirano, E., & Almstedt, A. E. (1996).
        Eulerian two-phase flow theory applied to fluidisation.
        International Journal of Multiphase Flow, 22, 21-66.
    \endverbatim

SourceFiles
    Ergun.C

\*---------------------------------------------------------------------------*/

#ifndef Ergun_H
#define Ergun_H

#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                            Class Ergun Declaration
\*---------------------------------------------------------------------------*/

class Ergun
:
    public dispersedDragModel
{
public:

    //- Runtime type information
    TypeName("Ergun");


    // Constructors

        //- Construct from a dictionary and an interface
        Ergun
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~Ergun();


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
