/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::functionObjects::pressure

Description
    Includes tools to manipulate the pressure into different forms.

    These currently include:
    - static pressure
        \f[
            p = \rho p_k
        \f]
    - total pressure
        \f[
            p_0 = p_{ref} + p + 0.5 \rho |U|^2
        \f]
    - static pressure coefficient
        \f[
            Cp = \frac{p - p_{\inf}}{0.5 \rho_{\inf} |U_{\inf}|^2}
        \f]
    - total pressure coefficient
        \f[
            Cp_0 = \frac{p_0 - p_{\inf}}{0.5 \rho_{\inf} |U_{\inf}|^2}
        \f]

    where
    \vartable
        \rho        | Density [kg/m^3]
        U           | Velocity [m/s]
        \rho_{\inf} | Freestream density [kg/m^3]
        p_{\inf}    | Freestream pressure [Pa]
        U_{\inf}    | Freestream velocity [m/s]
        p_k         | Kinematic pressure (p/rho) [m^2/s^2]
        p           | Pressure [Pa]
        p_0         | Total pressure [Pa]
        p_{ref}     | Reference pressure level [Pa]
        Cp          | Pressure coefficient
        Cp_0        | Total pressure coefficient
    \endvartable

    The function object will operate on both kinematic (\f$ p_k \f$) and static
    pressure (\f$ p \f$) fields, and the result is written as a
    volScalarField.

    The modes of operation are:
    \table
        Mode                        | calcTotal | calcCoeff
        Static pressure             | no        | no
        Total pressure              | yes       | no
        Pressure coefficient        | no        | yes
        Total pressure coefficient  | yes       | yes
    \endtable

    Example of function object specification to calculate pressure coefficient:
    \verbatim
    pressure1
    {
        type        pressure;
        libs        ("libfieldFunctionObjects.so");
        ...
        calcTotal   no;
        calcCoeff   yes;
    }
    \endverbatim

Usage
    \table
        Property     | Description                 | Required   | Default value
        type         | type name: pressure    | yes        |
        field        | Name of the pressure field  | no         | p
        U            | Name of the velocity field  | no         | U
        rho          | Name of the density field   | no         | rho
        result       | Name of the resulting field | no         | derived from p
        calcTotal    | Calculate total coefficient | yes        |
        pRef         | Reference pressure for total pressure | no | 0
        calcCoeff    | Calculate pressure coefficient | yes  |
        pInf         | Freestream pressure for coefficient calculation | no |
        UInf         | Freestream velocity for coefficient calculation | no |
        rhoInf       | Freestream density for coefficient calculation  | no |
    \endtable

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    pressure.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_pressure_H
#define functionObjects_pressure_H

#include "fieldExpression.H"
#include "volFieldsFwd.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class pressure Declaration
\*---------------------------------------------------------------------------*/

class pressure
:
    public fieldExpression
{
    // Private Data

        //- Name of velocity field, default is "U"
        word UName_;

        //- Name of density field, default is "rho"
        word rhoName_;


        // Total pressure calculation

            //- Flag to calculate total pressure
            bool calcTotal_;

            //- Reference pressure level
            scalar pRef_;


        // Pressure coefficient calculation

            //- Flag to calculate pressure coefficient
            bool calcCoeff_;

            //- Freestream pressure
            scalar pInf_;

            //- Freestream velocity
            vector UInf_;

            //- Freestream density
            scalar rhoInf_;


    // Private Member Functions

        //- Return the name of the derived pressure field
        word resultName() const;

        //- Multiply the static pressure p by rhoInf if necessary and return
        tmp<volScalarField> rhoScale(const volScalarField& p) const;

        //- Multiply the given field by rho or rhoInf as appropriate and return
        tmp<volScalarField> rhoScale
        (
            const volScalarField& p,
            const tmp<volScalarField>& tsf
        ) const;

        //- Return the reference pressure
        tmp<volScalarField> pRef(const tmp<volScalarField>& tp) const;

        //- Calculate and return the dynamic pressure
        tmp<volScalarField> pDyn
        (
            const volScalarField& p,
            const tmp<volScalarField>& tp
        ) const;

        //- Convert to coeff by applying the freestream dynamic pressure scaling
        tmp<volScalarField> coeff(const tmp<volScalarField>& tp) const;

        //- Calculate the derived pressure field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("pressure");


    // Constructors

        //- Construct from Time and dictionary
        pressure
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~pressure();


    // Member Functions

        //- Read the pressure data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
