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
    Foam::contactAngleModels::gravitational

Description
    Gravitational acceleration based contact angle model.

Usage:
    \table
        Property | Description                      | Required | Default value
        theta0   | Equilibrium contact angle        | yes      |
        gTheta   | Gravitational acceleration scale | yes      |
        thetaAdv | Limiting advancing contact angle | yes      |
        thetaRec | Limiting receding contact angle  | yes      |
    \endtable

    Example:
    \verbatim
    contactAngle
    {
        type           gravitational;
        gTheta         0.1;
        theta0         70;
        thetaAdv       100;
        thetaRec       50;
    }
    \endverbatim

See also
    Foam::contactAngleModel

SourceFiles
    gravitationalContactAngle.C

\*---------------------------------------------------------------------------*/

#ifndef gravitationalContactAngle_H
#define gravitationalContactAngle_H

#include "contactAngleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace contactAngleModels
{

/*---------------------------------------------------------------------------*\
                          Class gravitational Declaration
\*---------------------------------------------------------------------------*/

class gravitational
:
    public contactAngleModel
{
    // Private Data

        //- Equilibrium contact angle
        scalar theta0_;

        //- Gavitational acceleration scale
        scalar gTheta_;

        //- Limiting advancing contact angle
        scalar thetaAdv_;

        //- Limiting receding contact angle
        scalar thetaRec_;


public:

    //- Runtime type information
    TypeName("gravitational");


    // Constructors

        //- Construct from dictionary
        gravitational(const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<contactAngleModel> clone() const
        {
            return autoPtr<contactAngleModel>(new gravitational(*this));
        }


    //- Destructor
    virtual ~gravitational();


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
