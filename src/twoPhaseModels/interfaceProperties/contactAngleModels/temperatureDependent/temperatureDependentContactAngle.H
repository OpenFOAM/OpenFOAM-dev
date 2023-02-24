/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::contactAngleModels::temperatureDependent

Description
    Temperature-dependent contact angle model.

    The contact angle is evaluated from the specified Foam::Function1 for the
    temperature field looked-up from the mesh database the name of which
    may optionally be provided.

Usage
    \table
        Property     | Description               | Required    | Default value
        T            | Temperature field name    | no          | T
        contactAngle | contact angle function  | yes         |
    \endtable

    Example of the contact angle specification:
    \verbatim
        contactAngle
        {
            type                temperatureDependent;
            contactAngle        constant 0.07;
        }
    \endverbatim

See also
    Foam::contactAngleModel
    Foam::Function1

SourceFiles
    temperatureDependentContactAngle.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDependentContactAngle_H
#define temperatureDependentContactAngle_H

#include "contactAngleModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace contactAngleModels
{

/*---------------------------------------------------------------------------*\
                    Class temperatureDependent Declaration
\*---------------------------------------------------------------------------*/

class temperatureDependent
:
    public contactAngleModel
{
    // Private Data

        //- Name of temperature field, default = "T"
        word TName_;

        //- Equilibrium contact angle function
        autoPtr<Function1<scalar>> theta0_;


public:

    //- Runtime type information
    TypeName("temperatureDependent");


    // Constructors

        //- Construct from dictionary
        temperatureDependent(const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<contactAngleModel> clone() const
        {
            return autoPtr<contactAngleModel>(new temperatureDependent(*this));
        }


    //- Destructor
    virtual ~temperatureDependent();


    // Member Functions

        //- Cosine of the contact angle
        virtual tmp<scalarField> cosTheta
        (
            const fvPatchVectorField& Up,
            const vectorField& nHat
        ) const;

        //- Write data in dictionary format
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace contactAngleModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
