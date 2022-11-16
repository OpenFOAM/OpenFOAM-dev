/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::dispersedTurbulentDispersionModel

Description
    Model for turbulent dispersion between two phases where one phase can be
    considered dispersed in the other

SourceFiles
    dispersedTurbulentDispersionModel.C

\*---------------------------------------------------------------------------*/

#ifndef dispersedTurbulentDispersionModel_H
#define dispersedTurbulentDispersionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "turbulentDispersionModel.H"
#include "dispersedPhaseInterface.H"
#include "phaseCompressibleMomentumTransportModel.H"

namespace Foam
{
namespace turbulentDispersionModels
{

/*---------------------------------------------------------------------------*\
             Class dispersedTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class dispersedTurbulentDispersionModel
:
    public turbulentDispersionModel
{
protected:

    // Protected data

        //- Interface
        const dispersedPhaseInterface interface_;


public:

    // Constructors

        // Construct from a dictionary and an interface
        dispersedTurbulentDispersionModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    //- Destructor
    virtual ~dispersedTurbulentDispersionModel();


    // Member Functions

        //- Return a reference to the turbulence model for the continuous phase
        const phaseCompressible::momentumTransportModel&
            continuousTurbulence() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
