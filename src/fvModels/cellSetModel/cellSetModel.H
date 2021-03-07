/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Foam::fv::cellSetModel

Description
    Cell-set fvModels abstract base class. Provides a base set of controls
    regarding the location where the fvModel is applied.

Usage
    Example usage:
    \verbatim
    fvModel1
    {
        type            <fvModelType>

        // Apply everywhere
        selectionMode   all;

        // // Apply within a given cell set
        // selectionMode   cellSet;
        // cellSet         c0;

        // // Apply in cells containing a list of points
        // selectionMode   points;
        // points
        // (
        //     (2.25 0.5 0)
        //     (2.75 0.5 0)
        // );

        ...
    }
    \endverbatim

SourceFiles
    cellSetModel.C

\*---------------------------------------------------------------------------*/

#ifndef cellSetModel_H
#define cellSetModel_H

#include "fvModel.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class cellSetModel Declaration
\*---------------------------------------------------------------------------*/

class cellSetModel
:
    public fvModel
{
public:

    // Public data

        //- Enumeration for selection mode types
        enum class selectionModeType
        {
            points,
            cellSet,
            cellZone,
            all
        };

        //- Word list of selection mode type names
        static const NamedEnum<selectionModeType, 4>
            selectionModeTypeNames_;


private:

    // Private data

        //- Cell selection mode
        selectionModeType selectionMode_;

        //- Name of cell set for "cellSet" and "cellZone" selectionMode
        word cellSetName_;

        //- List of points for "points" selectionMode
        List<point> points_;

        //- Set of cells to apply source to
        mutable labelList cells_;

        //- Sum of cell volumes
        mutable scalar V_;


    // Private functions

        //- Non-virtual read
        void readCoeffs();

        //- Set the cell set based on the user input selection mode
        void setCellSet();


public:

    //- Runtime type information
    TypeName("cellSetModel");


    // Constructors

        //- Construct from components
        cellSetModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~cellSetModel();


    // Member Functions

        // Access

            //- Return const access to the cell selection mode
            inline const selectionModeType& selectionMode() const;

            //- Return const access to the name of cell set for "cellSet"
            //  selectionMode
            inline const word& cellSetName() const;

            //- Return const access to the total cell volume
            inline scalar V() const;

            //- Return const access to the cell set
            inline const labelList& cells() const;


        // Mesh changes

            //- Update for mesh changes
            virtual void updateMesh(const mapPolyMesh&);

            //- Update for mesh motion
            virtual bool movePoints();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellSetModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
