/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::fv::damping

Description
    Base fvModel for damping functions.

See also
    Foam::fv::isotropicDamping
    Foam::fv::verticalDamping

SourceFiles
    damping.C

\*---------------------------------------------------------------------------*/

#ifndef damping_H
#define damping_H

#include "fvModel.H"
#include "Function1.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class damping Declaration
\*---------------------------------------------------------------------------*/

class damping
:
    public fvModel
{
protected:

    // Protected Data

        //- Name of the velocity field
        word UName_;

        //- Damping coefficient [1/s]
        dimensionedScalar lambda_;

        //- The scaling function
        autoPtr<Function1<scalar>> scale_;

        //- Origins of the scaling coordinate
        vectorField origins_;

        //- Directions of increasing scaling coordinate
        vectorField directions_;


    // Protected Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Return the force coefficient
        tmp<volScalarField::Internal> forceCoeff() const;


public:

    //- Runtime type information
    TypeName("damping");


    // Constructors

        //- Construct from components
        damping
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~damping()
    {}


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
