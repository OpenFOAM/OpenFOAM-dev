/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::functionObjects::forcesBase

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::timeControl

SourceFiles
    forcesBase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forcesBase_H
#define functionObjects_forcesBase_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "volFields.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class forcesBase Declaration
\*---------------------------------------------------------------------------*/

class forcesBase
:
    public fvMeshFunctionObject,
    public logFiles
{

protected:

    // Protected data

        //- Enumeration for ensuring the right file is accessed
        enum class fileID
        {
            mainFile = 0,
            binsFile = 1
        };

        //- Pressure, viscous and porous force per bin
        List<Field<vector>> force_;

        //- Pressure, viscous and porous moment per bin
        List<Field<vector>> moment_;


        // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of density field (optional)
            word rhoName_;

            //- The name of the phase (optional)
            word phaseName_;

            //- Is the force density being supplied directly?
            Switch directForceDensity_;

            //- The name of the force density (fD) field
            word fDName_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;

            //- Flag to include porosity effects
            bool porosity_;


            // Bin information

                //- Number of bins
                label nBin_;

                //- Direction used to determine bin orientation
                vector binDir_;

                //- Distance between bin divisions
                scalar binDx_;

                //- Minimum bin bounds
                scalar binMin_;

                //- Bin positions along binDir
                List<point> binPoints_;

                //- Should bin data be cumulative?
                bool binCumulative_;


            //- Initialised flag
            bool initialised_;


    // Protected Member Functions

        using logFiles::file;

        Ostream& file(const fileID fid)
        {
            return logFiles::file(label(fid));
        }

        //- Create file names for forces and bins
        wordList createFileNames(const dictionary& dict) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);

        //- Initialise the fields
        void initialise();

        //- Return the effective surface stress
        tmp<surfaceVectorField> devTau() const;

        //- Dynamic viscosity field
        tmp<volScalarField> mu() const;

        //- Return rho if specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
        //  otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- Get the volume fraction field
        tmp<volScalarField> alpha() const;

        //- Get the volume fraction field on a patch
        tmp<scalarField> alpha(const label patchi) const;

        //- Accumulate bin data
        void applyBins
        (
            const vectorField& Md,
            const vectorField& fN,
            const vectorField& fT,
            const vectorField& fP,
            const vectorField& d
        );

        //- Calculate the forces and moments
        void calcForcesMoments(const vector& CofR);

        //- Return the current centre of the rigid body
        virtual vector CofR() const = 0;

        //- Write the constant centre of rotation value in the header
        //  Not written by default
        virtual void writeCoRValueHeader(Ostream& file);

        //- Write the time varying centre of rotation column header
        virtual void writeCoRHeader(Ostream& file);

        //- Write the time varying centre of rotation
        virtual void writeCofR(Ostream& file);

        //- Helper function to write force data
        void writeForces();

        //- Helper function to write bin data
        void writeBins();


public:

    //- Runtime type information
    TypeName("forcesBase");


    // Constructors

        //- Construct from Time and dictionary
        forcesBase
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from objectRegistry and dictionary
        forcesBase
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        forcesBase(const forcesBase&) = delete;


    //- Destructor
    virtual ~forcesBase();


    // Member Functions

        //- Read the forces data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Calculate the forces and moments
        virtual void calcForcesMoments();

        //- Return the total force
        virtual vector forceEff() const;

        //- Return the total moment
        virtual vector momentEff() const;

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the forces
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const forcesBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
