/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2020 OpenFOAM Foundation
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
    Foam::phaseTransferModel

Description

SourceFiles
    phaseTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseTransferModel_H
#define phaseTransferModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "HashPtrTable.H"
#include "hashedWordList.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                     Class phaseTransferModel Declaration
\*---------------------------------------------------------------------------*/

class phaseTransferModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("phaseTransferModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseTransferModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static Data Members

        //- Mass transfer rate dimensions
        static const dimensionSet dimDmdt;

        //- Empty species list
        static const hashedWordList noSpecies_;


    // Constructors

        //- Construct from a dictionary and a phase pair
        phaseTransferModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~phaseTransferModel();


    // Selectors

        static autoPtr<phaseTransferModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Whether or not there is a mass transfer rate for the mixture
        virtual bool mixture() const;

        //- The mass transfer rate for the mixture
        virtual tmp<volScalarField> dmdtf() const;

        //- The list of individual species that are transferred
        virtual const hashedWordList& species() const;

        //- The mass transfer rate for individual species
        virtual HashPtrTable<volScalarField> dmidtf() const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
