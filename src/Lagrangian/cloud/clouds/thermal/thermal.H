/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::clouds::thermal

Description
    Base class for clouds with thermodynamic modelling

SourceFiles
    thermal.C

\*---------------------------------------------------------------------------*/

#ifndef thermal_H
#define thermal_H

#include "massive.H"
#include "carried.H"
#include "basicLagrangianThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                             Class thermal Declaration
\*---------------------------------------------------------------------------*/

class thermal
:
    private autoPtr<basicLagrangianThermo>,
    public massive
{
private:

    // Private Data

        //- Reference to the cloud
        const cloud& cloud_;


protected:

    // Protected Constructors

        //- Construct from a reference to the cloud and its base classes and a
        //  thermo model pointer
        thermal(const cloud&, const shaped&, basicLagrangianThermo*);


public:

    // Public Static Data

        //- Run-time type information
        TypeName("thermal");


    // Constructors

        //- Construct from a reference to the cloud and its base classes
        thermal(const cloud&, const shaped&, const carried&);


    // Public Data

        //- Temperature
        const LagrangianScalarDynamicField& T;

        //- Internal energy
        CloudStateFieldRef<scalar> e;

        //- Heat capacity at constant volume
        const LagrangianScalarDynamicField& Cv;

        //- Thermal conductivity
        const LagrangianScalarDynamicField& kappa;


    //- Destructor
    virtual ~thermal();


    // Member Functions

        //- Determine whether the thermo model is one of the given types
        template<class Thermo, class ... Thermos>
        bool isThermo() const;

        //- Generate an error if the thermo model is not one of the given types
        template<class Thermo, class ... Thermos>
        void assertThermo() const;

        //- Generate an error if the thermo model is not one of the given
        //  types. A model is provided to give the error message some context.
        template<class Thermo, class ... Thermos>
        void assertThermo(const LagrangianModel& model) const;

        //- Get a reference to the thermodynamic model of a given type
        template<class Thermo = basicLagrangianThermo, class ... Args>
        const Thermo& thermo(const Args& ... args) const;

        //- Non-const override of the above
        template<class Thermo = basicLagrangianThermo, class ... Args>
        Thermo& thermo(const Args& ... args);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thermalTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
