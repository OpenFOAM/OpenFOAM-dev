/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::combustionModels::zoneCombustion

Description
    Zone-filtered combustion model.

    Enable the reactions within the specified list of cell-zones and set
    to zero elsewhere.

SourceFiles
    zoneCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef zoneCombustion_H
#define zoneCombustion_H

#include "CombustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class zoneCombustion Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo>
class zoneCombustion
:
    public CombustionModel<ReactionThermo>
{
    // Private data

        //- The combustion model to be zone-filtered
        autoPtr<CombustionModel<ReactionThermo>> combustionModelPtr_;

        //- List of zone names in which the reactions are active
        wordList zoneNames_;


    // Private Member Functions

        //- Filter the reaction-rate matrix on the cellZones
        tmp<fvScalarMatrix> filter(const tmp<fvScalarMatrix>& tR) const;

        //- Filter the given field on the cellZones
        tmp<volScalarField> filter(const tmp<volScalarField>& tS) const;

        //- Disallow copy construct
        zoneCombustion(const zoneCombustion&);

        //- Disallow default bitwise assignment
        void operator=(const zoneCombustion&);


public:

    //- Runtime type information
    TypeName("zoneCombustion");


    // Constructors

        //- Construct from components
        zoneCombustion
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~zoneCombustion();


    // Member Functions

        //- Return access to the thermo package
        virtual ReactionThermo& thermo();

        //- Return const access to the thermo package
        virtual const ReactionThermo& thermo() const;

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "zoneCombustion.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
