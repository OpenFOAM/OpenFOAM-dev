/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
    Foam::saturationPressureModel

Description
    Model to describe the dependence of saturation pressure on temperature, and
    vice versa

SourceFiles
    saturationPressureModel.C
    saturationPressureModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef saturationPressureModel_H
#define saturationPressureModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DEFINE_PSAT(FieldType, Modifier)                                       \
                                                                               \
    virtual tmp<FieldType> pSat(const FieldType& T) const Modifier;            \
                                                                               \
    virtual tmp<FieldType> pSatPrime(const FieldType& T) const Modifier;       \
                                                                               \
    virtual tmp<FieldType> lnPSat(const FieldType& T) const Modifier;


#define IMPLEMENT_PSAT(ModelType, FieldType)                                   \
                                                                               \
    Foam::tmp<Foam::FieldType>                                                 \
    Foam::saturationModels::ModelType::pSat(const FieldType& T) const          \
    {                                                                          \
        return pSat<FieldType>(T);                                             \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::FieldType>                                                 \
    Foam::saturationModels::ModelType::pSatPrime(const FieldType& T) const     \
    {                                                                          \
        return pSatPrime<FieldType>(T);                                        \
    }                                                                          \
                                                                               \
    Foam::tmp<Foam::FieldType>                                                 \
    Foam::saturationModels::ModelType::lnPSat(const FieldType& T) const        \
    {                                                                          \
        return lnPSat<FieldType>(T);                                           \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class saturationPressureModel Declaration
\*---------------------------------------------------------------------------*/

class saturationPressureModel
{
public:

    //- Runtime type information
    TypeName("saturationPressureModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        saturationPressureModel,
        dictionary,
        (const dictionary& dict),
        (dict)
    );


    // Constructors

        //- Default construct
        saturationPressureModel();


    // Selectors

        //- Select with dictionary
        static autoPtr<saturationPressureModel> New
        (
            const dictionary& dict
        );

        //- Select with name within a dictionary
        static autoPtr<saturationPressureModel> New
        (
            const word& name,
            const dictionary& dict
        );


    //- Destructor
    virtual ~saturationPressureModel();


    // Member Functions

        //- Saturation pressure for volScalarField::Internal
        DEFINE_PSAT(volScalarField::Internal, = 0);

        //- Saturation pressure for volScalarField
        DEFINE_PSAT(volScalarField, = 0);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
