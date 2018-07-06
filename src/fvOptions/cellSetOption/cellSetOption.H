/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fv::cellSetOption

Description
    Cell-set options abtract base class.  Provides a base set of controls,
    e.g.:
    \verbatim
        type            scalarExplicitSource    // Source type
        active          on;                     // on/off switch

        timeStart       0.0;        // Start time
        duration        1000.0;     // Duration
        selectionMode   cellSet;    // cellSet, points, cellZone
        .
        .
        .
    \endverbatim

Note
    Source/sink options are to be added to the equation R.H.S.

SourceFiles
    cellSetOption.C
    cellSetOptionIO.C

\*---------------------------------------------------------------------------*/

#ifndef cellSetOption_H
#define cellSetOption_H

#include "fvOption.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class cellSetOption Declaration
\*---------------------------------------------------------------------------*/

class cellSetOption
:
    public option
{
public:

    // Public data

        //- Enumeration for selection mode types
        enum selectionModeType
        {
            smPoints,
            smCellSet,
            smCellZone,
            smAll
        };

        //- Word list of selection mode type names
        static const NamedEnum<selectionModeType, 4>
            selectionModeTypeNames_;


protected:

    // Protected data

        //- Time start
        scalar timeStart_;

        //- Duration
        scalar duration_;

        //- Cell selection mode
        selectionModeType selectionMode_;

        //- Name of cell set for "cellSet" and "cellZone" selectionMode
        word cellSetName_;

        //- List of points for "points" selectionMode
        List<point> points_;

        //- Set of cells to apply source to
        labelList cells_;

        //- Sum of cell volumes
        scalar V_;


    // Protected functions

        //- Set the cellSet or points selection
        void setSelection(const dictionary& dict);

        //- Set the cell set based on the user input selection mode
        void setCellSet();


public:

    //- Runtime type information
    TypeName("cellSetOption");


    // Constructors

        //- Construct from components
        cellSetOption
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~cellSetOption();


    // Member Functions

        // Access

            //- Return const access to the time start
            inline scalar timeStart() const;

            //- Return const access to the duration
            inline scalar duration() const;

            //- Return true if within time limits
            inline bool inTimeLimits(const scalar time) const;

            //- Return const access to the cell selection mode
            inline const selectionModeType& selectionMode() const;

            //- Return const access to the name of cell set for "cellSet"
            //  selectionMode
            inline const word& cellSetName() const;

            //- Return const access to the total cell volume
            inline scalar V() const;

            //- Return const access to the cell set
            inline const labelList& cells() const;


        // Edit

            //- Return access to the time start
            inline scalar& timeStart();

            //- Return access to the duration
            inline scalar& duration();


        // Checks

            //- Is the source active?
            virtual bool isActive();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellSetOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
