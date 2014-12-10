/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::MRFZone

Description
    MRF zone definition based on cell zone and parameters
    obtained from a control dictionary constructed from the given stream.

    The rotation of the MRF region is defined by an origin and axis of
    rotation and an angular speed.

SourceFiles
    MRFZone.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZone_H
#define MRFZone_H

#include "dictionary.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "fvMatricesFwd.H"
#include "mapPolyMesh.H"
#include "DataEntry.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class MRFZone Declaration
\*---------------------------------------------------------------------------*/

class MRFZone
{
    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Name of the MRF region
        const word name_;

        //- Coefficients dictionary
        dictionary coeffs_;

        //- MRF region active flag
        bool active_;

        //- Name of cell zone
        word cellZoneName_;

        //- Cell zone ID
        label cellZoneID_;

        const wordList excludedPatchNames_;
        labelList excludedPatchLabels_;

        //- Internal faces that are part of MRF
        labelList internalFaces_;

        //- Outside faces (per patch) that move with the MRF
        labelListList includedFaces_;

        //- Excluded faces (per patch) that do not move with the MRF
        labelListList excludedFaces_;

        //- Origin of the axis
        const vector origin_;

        //- Axis vector
        vector axis_;

        //- Angular velocty (rad/sec)
        autoPtr<DataEntry<scalar> > omega_;


    // Private Member Functions

        //- Divide faces in frame according to patch
        void setMRFFaces();

        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType>
        void makeRelativeRhoFlux
        (
            const RhoFieldType& rho,
            surfaceScalarField& phi
        ) const;

        //- Make the given absolute mass/vol flux relative within the MRF region
        template<class RhoFieldType>
        void makeRelativeRhoFlux
        (
            const RhoFieldType& rho,
            FieldField<fvsPatchField, scalar>& phi
        ) const;

        //- Make the given relative mass/vol flux absolute within the MRF region
        template<class RhoFieldType>
        void makeAbsoluteRhoFlux
        (
            const RhoFieldType& rho,
            surfaceScalarField& phi
        ) const;

        //- Disallow default bitwise copy construct
        MRFZone(const MRFZone&);

        //- Disallow default bitwise assignment
        void operator=(const MRFZone&);


public:

    // Declare name of the class and its debug switch
    ClassName("MRFZone");


    // Constructors

        //- Construct from fvMesh
        MRFZone
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& cellZoneName = word::null
        );

        //- Return clone
        autoPtr<MRFZone> clone() const
        {
            notImplemented("autoPtr<MRFZone> clone() const");
            return autoPtr<MRFZone>(NULL);
        }


    // Member Functions

        // Access

            //- Return const access to the MRF region name
            inline const word& name() const;

            //- Return const access to the MRF active flag
            inline bool active() const;

            //- Return the current Omega vector
            vector Omega() const;


        // Evaluation

            //- Update the mesh corresponding to given map
            void updateMesh(const mapPolyMesh& mpm)
            {
                // Only updates face addressing
                setMRFFaces();
            }

            //- Add the Coriolis force contribution to the acceleration field
            void addCoriolis
            (
                const volVectorField& U,
                volVectorField& ddtU
            ) const;

            //- Add the Coriolis force contribution to the momentum equation
            //  Adds to the lhs of the equation; optionally add to rhs
            void addCoriolis
            (
                fvVectorMatrix& UEqn,
                const bool rhs = false
            ) const;

            //- Add the Coriolis force contribution to the momentum equation
            //  Adds to the lhs of the equation; optionally add to rhs
            void addCoriolis
            (
                const volScalarField& rho,
                fvVectorMatrix& UEqn,
                const bool rhs = false
            ) const;

            //- Make the given absolute velocity relative within the MRF region
            void makeRelative(volVectorField& U) const;

            //- Make the given absolute flux relative within the MRF region
            void makeRelative(surfaceScalarField& phi) const;

            //- Make the given absolute boundary flux relative
            //  within the MRF region
            void makeRelative(FieldField<fvsPatchField, scalar>& phi) const;

            //- Make the given absolute mass-flux relative within the MRF region
            void makeRelative
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;

            //- Make the given relative velocity absolute within the MRF region
            void makeAbsolute(volVectorField& U) const;

            //- Make the given relative flux absolute within the MRF region
            void makeAbsolute(surfaceScalarField& phi) const;

            //- Make the given relative mass-flux absolute within the MRF region
            void makeAbsolute
            (
                const surfaceScalarField& rho,
                surfaceScalarField& phi
            ) const;

            //- Correct the boundary velocity for the roation of the MRF region
            void correctBoundaryVelocity(volVectorField& U) const;


    // I-O

        //- Write
        void writeData(Ostream& os) const;

        //- Read MRF dictionary
        bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MRFZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MRFZoneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
