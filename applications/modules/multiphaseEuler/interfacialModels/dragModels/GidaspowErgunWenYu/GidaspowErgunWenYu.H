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
    Foam::dragModels::GidaspowErgunWenYu

Description
    Gidaspow, Ergun, Wen and Yu drag model

    Reference:
    \verbatim
        Gidaspow, D. (1994).
        Multiphase flow and fluidisation: continuum and kinetic theory
        descriptions.
        Academic press, New York.
    \endverbatim

SourceFiles
    GidaspowErgunWenYu.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowErgunWenYu_H
#define GidaspowErgunWenYu_H

#include "dispersedDragModel.H"
#include "Ergun.H"
#include "WenYu.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                     Class GidaspowErgunWenYu Declaration
\*---------------------------------------------------------------------------*/

class GidaspowErgunWenYu
:
    public dispersedDragModel
{
    // Private Data

        //- Ergun drag model
        Ergun Ergun_;

        //- Wen Yu drag model
        WenYu WenYu_;


public:

    //- Runtime type information
    TypeName("GidaspowErgunWenYu");


    // Constructors

        //- Construct from a dictionary and an interface
        GidaspowErgunWenYu
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool registerObject
        );


    //- Destructor
    virtual ~GidaspowErgunWenYu();


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
