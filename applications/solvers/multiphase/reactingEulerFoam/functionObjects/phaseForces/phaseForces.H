/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::functionObjects::phaseForces

Group
    grpForcesFunctionObjects

Description
    This function object calculates and outputs the blended interfacial forces
    acting on a given phase (drag, virtual mass, lift, wall-lubrication and
    turbulent dispersion). It works only in combination with the
    reactingEulerFoam solvers. For a simulation involving more than two phases,
    the accumulated force is calculated by looping over all phasePairs involving
    that phase. The fields are stored in the database so that they can be
    processed further, e.g. with the fieldAveraging functionObject.

    Example of function object specification:
    \verbatim
    phaseForces.water
    {
        type            phaseForces;
        libs ("libreactingEulerFoamFunctionObjects.so");
        writeControl    outputTime;
        writeInterval   1;
        log             false;
        ...
        phaseName       water;
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type name: phaseForces  | yes         |
        phaseName    | Name of evaluated phase | yes         |
    \endtable

See also
    Foam::BlendedInterfacialModel
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    phaseForces.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_phaseForces_H
#define functionObjects_phaseForces_H

#include "fvMeshFunctionObject.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                         Class phaseForces Declaration
\*---------------------------------------------------------------------------*/

class phaseForces
:
    public fvMeshFunctionObject
{
    // Private data

        //- Phase for which forces are evaluated
        const phaseModel& phase_;

        //- Constant access to phaseSystem
        const phaseSystem& fluid_;

        //- Drag force
        volVectorField dragForce_;

        //- Virtual mass force
        volVectorField virtualMassForce_;

        //- Lift force
        volVectorField liftForce_;

        //- Wall lubrication force
        volVectorField wallLubricationForce_;

        //- Turbulent dispersion force
        volVectorField turbulentDispersionForce_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        phaseForces(const phaseForces&);

        //- Disallow default bitwise assignment
        void operator=(const phaseForces&);

        //- Evaluate and return the drag force
        tmp<volVectorField> dragForce() const;

        //- Evaluate and return the virtual mass force
        tmp<volVectorField> virtualMassForce() const;

        //- Evaluate and return the lift force
        tmp<volVectorField> liftForce() const;

        //- Evaluate and return the wall lubrication force
        tmp<volVectorField> wallLubricationForce() const;

        //- Evaluate and return the turbulent dispersion force
        tmp<volVectorField> turbulentDispersionForce() const;


public:

    //- Runtime type information
    TypeName("phaseForces");


    // Constructors

        //- Construct from Time and dictionary
        phaseForces
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~phaseForces();


    // Member Functions

        //- Read the input data
        virtual bool read(const dictionary& dict);

        //- Calculate the force fields
        virtual bool execute();

        //- Write the force fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
