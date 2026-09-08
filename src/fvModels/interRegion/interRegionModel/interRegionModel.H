/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::fv::interRegionModel

Description
    Base class for inter-region exchange.

\*---------------------------------------------------------------------------*/

#ifndef interRegionModel_H
#define interRegionModel_H

#include "fvModel.H"
#include "volFields.H"
#include "cellsToCells.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class interRegionModel Declaration
\*---------------------------------------------------------------------------*/

class interRegionModel
:
    public fvModel
{
    // Private data

        //- Is this the owner region? +1 is yes, 0 is no, -1 is either "yet to
        //  be determined" or indicates that this model has no neighbouring
        //  model (e.g., inter-region porosity)
        mutable label owner_;

        //- Name of the neighbour region to map
        word nbrRegionName_;

        //- Name of the corresponding model in the neighbour region
        word nbrModelName_;

        //- Cached dictionary for delayed construction of nbrModelName
        dictionary nbrModelNameDict_;

        //- Interpolation method
        word interpolationMethod_;

        //- Mesh to mesh interpolation object
        mutable autoPtr<cellsToCells> interpolationPtr_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Get the neighbour interRegionModel, or null if there is none
        const interRegionModel& nbrModelOrNull() const;

        //- Access the interpolation engine
        const cellsToCells& interpolation() const;


protected:

    // Protected Member Functions

        //- Get the neighbour interRegionModel
        template<class InterRegionModelType = interRegionModel>
        const InterRegionModelType& nbrModel() const;


public:

    //- Runtime type information
    TypeName("interRegionModel");


    // Constructors

        //- Construct from dictionary
        interRegionModel
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~interRegionModel();


    // Member Functions

        // Access

            //- Return whether the owner region
            inline bool owner() const;

            //- Return whether the neighbour region
            inline bool neighbour() const;

            //- Return const access to the neighbour region name
            inline const word& nbrRegionName() const;

            //- Return const access to the neighbour region name
            inline const word& nbrModelName() const;

            //- Return const access to the neighbour mesh
            inline const fvMesh& nbrMesh() const;


        // Interpolation

            //- Interpolate field
            template<class Type>
            tmp<Field<Type>> interpolate
            (
                const Field<Type>& field
            ) const;

            //- Interpolate field
            template<class Type>
            void interpolate
            (
                const Field<Type>& field,
                Field<Type>& result
            ) const;


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interRegionModelI.H"

#ifdef NoRepository
    #include "interRegionModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
