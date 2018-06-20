/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::functionObjects::CourantNo

Description
    Calculates and outputs the Courant number as a volScalarField.  The field is
    stored on the mesh database so that it can be retrieved and used for other
    applications.

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    CourantNo.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_CourantNo_H
#define functionObjects_CourantNo_H

#include "fieldExpression.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class CourantNo Declaration
\*---------------------------------------------------------------------------*/

class CourantNo
:
    public fieldExpression
{
    // Private data

        //- Name of flux field, default is "phi"
        word phiName_;

        //- Name of density field (optional)
        word rhoName_;


    // Private Member Functions

        //- Divide the Courant number by rho if required
        tmp<volScalarField::Internal> byRho
        (
            const tmp<volScalarField::Internal>& Co
        ) const;

        //- Calculate the Courant number field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("CourantNo");


    // Constructors

        //- Construct from Time and dictionary
        CourantNo
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    //- Destructor
    virtual ~CourantNo();


    // Member Functions

        //- Read the CourantNo data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
