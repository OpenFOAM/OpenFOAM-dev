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
    Foam::contactAngleModels::constant

Description
    Uniform constant contact angle model.

Usage
    Example of the contact angle specification:
    \verbatim
        contactAngle
        {
            type                constant;
            theta0              70;
        }
    \endverbatim

See also
    Foam::contactAngleModel

SourceFiles
    constantContactAngle.C

\*---------------------------------------------------------------------------*/

#ifndef constantContactAngle_H
#define constantContactAngle_H

#include "contactAngleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace contactAngleModels
{

/*---------------------------------------------------------------------------*\
                          Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public contactAngleModel
{
    // Private Data

        //- Equilibrium contact angle
        scalar theta0_;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from dictionary
        constant(const dictionary& dict);

        //- Construct and return a clone
        virtual autoPtr<contactAngleModel> clone() const
        {
            return autoPtr<contactAngleModel>(new constant(*this));
        }


    //- Destructor
    virtual ~constant();


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
