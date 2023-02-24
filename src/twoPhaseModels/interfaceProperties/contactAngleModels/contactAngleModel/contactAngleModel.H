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
    Foam::contactAngleModel

Description
    Abstract base-class for contact-angle models
    which return the cosine contact angle field.

Usage
    Example of the surface tension specification:
    \verbatim
        contactAngle
        {
            type                <contact angle model type>;
            <coefficient name>  <coefficient value>;
            .
            .
            .
        }
    \endverbatim

SourceFiles
    contactAngleModel.C
    contactAngleModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef contactAngleModel_H
#define contactAngleModel_H

#include "fvPatchFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class contactAngleModel Declaration
\*---------------------------------------------------------------------------*/

class contactAngleModel
{
protected:

    // Protected member functions

        static const dictionary& contactAngleDict(const dictionary& dict)
        {
            return dict.subDict("contactAngle");
        }


public:

    //- Runtime type information
    TypeName("contactAngleModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            contactAngleModel,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Constructors

        // Construct
        contactAngleModel();

        //- Construct and return a clone
        virtual autoPtr<contactAngleModel> clone() const = 0;


    //- Destructor
    virtual ~contactAngleModel();


    // Selectors

        static autoPtr<contactAngleModel> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Cosine of the contact angle
        virtual tmp<scalarField> cosTheta
        (
            const fvPatchVectorField& Up,
            const vectorField& nHat
        ) const = 0;

        //- Write data in dictionary format
        virtual void write(Ostream& os) const = 0;
};


void writeEntry(Ostream& os, const contactAngleModel& ca);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
