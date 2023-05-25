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
    Foam::dispersedLiftModel

Description
    Model for lift between two phases where one phase can be considered
    dispersed in the other and the lift therefore characterised by a lift
    coefficient

SourceFiles
    dispersedLiftModel.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedLiftModel_H
#define dispersedLiftModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "liftModel.H"
#include "dispersedPhaseInterface.H"

namespace Foam
{
namespace liftModels
{

/*---------------------------------------------------------------------------*\
                          Class dispersedLiftModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedLiftModel
:
    public liftModel
{
protected:

    // Protected data

        //- Interface
        const dispersedPhaseInterface interface_;


public:

    // Constructors

        // Construct from a dictionary and an interface
        dispersedLiftModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~dispersedLiftModel();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const = 0;

        //- Return phase-intensive lift force
        virtual tmp<volVectorField> Fi() const;

        //- Return lift force
        virtual tmp<volVectorField> F() const;

        //- Return face lift force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
