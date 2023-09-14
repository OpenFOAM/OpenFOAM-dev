/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::dispersedPhaseInterface

Description
    Class to represent a interface between phases where one phase is considered
    dispersed within the other.

SourceFiles
    dispersedPhaseInterface.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedPhaseInterface_H
#define dispersedPhaseInterface_H

#include "phaseInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class dispersedPhaseInterface Declaration
\*---------------------------------------------------------------------------*/

class dispersedPhaseInterface
:
    virtual public phaseInterface
{
private:

    // Private Data

        //- Dispersed phase
        const phaseModel& dispersed_;


public:

    //- Runtime type information
    TypeName("dispersedPhaseInterface");


    // Constructors

        //- Construct from phases
        dispersedPhaseInterface
        (
            const phaseModel& dispersed,
            const phaseModel& continuous
        );

        //- Construct from fluid and name
        dispersedPhaseInterface
        (
            const phaseSystem& fluid,
            const word& name
        );


    //- Destructor
    virtual ~dispersedPhaseInterface();


    // Static Member Functions

        //- Return the separator that delimits this interface's name
        static word separator()
        {
            return "dispersedIn";
        }


    // Member Functions

        //- Name
        virtual word name() const;


        // Access

            //- Dispersed phase
            const phaseModel& dispersed() const;

            //- Continuous phase
            const phaseModel& continuous() const;


        // Properties

            //- Relative velocity
            tmp<volVectorField> Ur() const;

            //- Reynolds number
            tmp<volScalarField> Re() const;

            //- Prandtl number
            tmp<volScalarField> Pr() const;

            //- Eotvos number
            tmp<volScalarField> Eo() const;

            //- Eotvos number for given diameter
            tmp<volScalarField> Eo(const volScalarField& d) const;

            //- Morton Number
            tmp<volScalarField> Mo() const;

            //- Takahashi Number
            tmp<volScalarField> Ta() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
