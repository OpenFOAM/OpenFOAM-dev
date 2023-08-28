/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class wordListAndType
{
    public:

        wordList wl;
        Type t;

        wordListAndType()
        {}

        wordListAndType(Istream& is)
        :
            wl(is),
            t(is)
        {}
};

template<class Type>
inline Istream& operator>>(Istream& is, wordListAndType<Type>& wlat)
{
    return is >> wlat.wl >> wlat.t;
}

template<class Type>
inline Ostream& operator<<(Ostream& os, const wordListAndType<Type>& wlat)
{
    return os << wlat.wl << wlat.t;
}

template<class Type>
inline bool operator==
(
    const wordListAndType<Type>& a,
    const wordListAndType<Type>& b
)
{
    return a.wl == b.wl && a.t == b.t;
}

template<class Type>
inline bool operator!=
(
    const wordListAndType<Type>& a,
    const wordListAndType<Type>& b
)
{
    return !(a == b);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::phaseSystem::mesh() const
{
    return mesh_;
}


inline const Foam::pimpleNoLoopControl& Foam::phaseSystem::pimple() const
{
    return pimple_;
}


inline const Foam::phaseSystem::phaseModelList&
Foam::phaseSystem::phases() const
{
    return phaseModels_;
}


inline Foam::phaseSystem::phaseModelList&
Foam::phaseSystem::phases()
{
    return phaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::movingPhases() const
{
    return movingPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::movingPhases()
{
    return movingPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::stationaryPhases() const
{
    return stationaryPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::stationaryPhases()
{
    return stationaryPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::anisothermalPhases() const
{
    return anisothermalPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::anisothermalPhases()
{
    return anisothermalPhaseModels_;
}


inline const Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::multicomponentPhases() const
{
    return multicomponentPhaseModels_;
}


inline Foam::phaseSystem::phaseModelPartialList&
Foam::phaseSystem::multicomponentPhases()
{
    return multicomponentPhaseModels_;
}


inline const Foam::phaseModel& Foam::phaseSystem::otherPhase
(
    const phaseModel& phase
) const
{
    if (phaseModels_.size() != 2)
    {
        FatalErrorInFunction
            << "Call from a two-phase model in a multi-phase system."
            << exit(FatalError);
    }

    if (&phase == &phaseModels_[0])
    {
        return phaseModels_[1];
    }
    else
    {
        return phaseModels_[0];
    }
}


inline const Foam::surfaceScalarField& Foam::phaseSystem::phi() const
{
    return phi_;
}


inline Foam::surfaceScalarField& Foam::phaseSystem::phi()
{
    return phi_;
}


inline const Foam::volScalarField& Foam::phaseSystem::dpdt() const
{
    return dpdt_;
}


inline Foam::volScalarField& Foam::phaseSystem::dpdt()
{
    return dpdt_;
}


inline const Foam::IOMRFZoneList& Foam::phaseSystem::MRF() const
{
    return MRF_;
}


inline Foam::fvModels& Foam::phaseSystem::fvModels(fvMesh& mesh)
{
    return Foam::fvModels::New(mesh);
}


inline const Foam::fvModels& Foam::phaseSystem::fvModels() const
{
    return Foam::fvModels::New(mesh_);
}


inline Foam::fvConstraints& Foam::phaseSystem::fvConstraints(fvMesh& mesh)
{
    return Foam::fvConstraints::New(mesh);
}


inline const Foam::fvConstraints& Foam::phaseSystem::fvConstraints() const
{
    return Foam::fvConstraints::New(mesh_);
}


inline const Foam::dimensionedScalar& Foam::phaseSystem::deltaN() const
{
    return deltaN_;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
inline void addField
(
    const label phasei,
    const word& name,
    tmp<GeoField> field,
    PtrList<GeoField>& fieldList
)
{
    if (fieldList.set(phasei))
    {
        fieldList[phasei] += field;
    }
    else
    {
        fieldList.set
        (
            phasei,
            new GeoField
            (
                name, // IOobject::groupName(name, group.name()),
                field
            )
        );
    }
}


template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    tmp<GeoField> field,
    PtrList<GeoField>& fieldList
)
{
    if (fieldList.set(group.index()))
    {
        fieldList[group.index()] += field;
    }
    else
    {
        fieldList.set
        (
            group.index(),
            new GeoField
            (
                IOobject::groupName(name, group.name()),
                field
            )
        );
    }
}


template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    const GeoField& field,
    PtrList<GeoField>& fieldList
)
{
    addField(group, name, tmp<GeoField>(field), fieldList);
}


template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    tmp<GeoField> field,
    HashPtrTable<GeoField>& fieldTable
)
{
    if (fieldTable.found(group.name()))
    {
        *fieldTable[group.name()] += field;
    }
    else
    {
        fieldTable.set
        (
            group.name(),
            new GeoField
            (
                IOobject::groupName(name, group.name()),
                field
            )
        );
    }
}


template<class GeoField, class Group>
inline void addField
(
    const Group& group,
    const word& name,
    const GeoField& field,
    HashPtrTable<GeoField>& fieldTable
)
{
    addField(group, name, tmp<GeoField>(field), fieldTable);
}


template<class Type, template<class> class PatchField, class GeoMesh>
PtrList<GeometricField<Type, PatchField, GeoMesh>> operator&
(
    const PtrList<PtrList<GeometricField<scalar, PatchField, GeoMesh>>>& As,
    const PtrList<GeometricField<Type, PatchField, GeoMesh>>& fs
)
{
    PtrList<GeometricField<Type, PatchField, GeoMesh >> Afs(fs.size());

    forAll(Afs, i)
    {
        Afs.set(i, As[i][i]*fs[i]);

        forAll(Afs, j)
        {
            if (j != i)
            {
                Afs[i] += As[i][j]*fs[j];
            }
        }
    }

    return Afs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
