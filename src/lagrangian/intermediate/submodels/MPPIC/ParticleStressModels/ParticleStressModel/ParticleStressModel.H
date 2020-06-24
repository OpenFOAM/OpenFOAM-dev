/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::ParticleStressModel

Description
    Base class for inter-particle stress models.

SourceFiles
    exponential.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleStressModel_H
#define ParticleStressModel_H

#include "Field.H"
#include "FieldField.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ParticleStressModel Declaration
\*---------------------------------------------------------------------------*/

class ParticleStressModel
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ParticleStressModel&) = delete;


protected:

    // Protected data

        //- Close pack volume fraction
        scalar alphaPacked_;


public:

    //- Runtime type information
    TypeName("particleStressModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ParticleStressModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    //- Constructors

        //- Construct from components
        ParticleStressModel(const dictionary& dict);

        //- Construct a copy
        ParticleStressModel(const ParticleStressModel& sm);

        //- Construct and return a clone
        virtual autoPtr<ParticleStressModel> clone() const = 0;


    //- Selector
    static autoPtr<ParticleStressModel> New
    (
        const dictionary& dict
    );


    //- Destructor
    virtual ~ParticleStressModel();


    //- Member Functions

        //- Access max volume fraction
        scalar alphaPacked() const;

        //- Collision stress
        virtual tmp<Field<scalar>> tau
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const = 0;

        //- Collision stress derivative w.r.t. the volume fraction
        virtual tmp<Field<scalar>> dTaudTheta
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const = 0;

        //- Collision stress using FieldFields
        tmp<FieldField<Field, scalar>> tau
        (
            const FieldField<Field, scalar>& alpha,
            const FieldField<Field, scalar>& rho,
            const FieldField<Field, scalar>& uRms
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
