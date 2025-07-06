/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

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

#include "basic2Thermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
/*
template<class Thermo, class Table>
typename Table::iterator Foam::basic2Thermo::lookupThermo
(
    const dictionary& thermoDict,
    Table* tablePtr //in v2412, this is a reference, not a pointer.
)
{
    word thermoTypeName;

    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        const int nCmpt = 7;
        const char* cmptNames[nCmpt] =
        {
            "type",
            "mixture",
            "transport",
            "thermo",
            "equationOfState",
            "specie",
            "energy"
        };

        // Construct the name of the thermo package from the components
        thermoTypeName =
            word(thermoTypeDict.lookup("type")) + '<'
          + word(thermoTypeDict.lookup("mixture")) + '<'
          + word(thermoTypeDict.lookup("transport")) + '<'
          + word(thermoTypeDict.lookup("thermo")) + '<'
          + word(thermoTypeDict.lookup("equationOfState")) + '<'
          + word(thermoTypeDict.lookup("specie")) + ">>,"
          + word(thermoTypeDict.lookup("energy")) + ">>>";

        // Lookup the thermo package
        //typename Table::iterator cstrIter = tablePtr->find(thermoTypeName);
        auto ctorIter = tablePtr->cfind(thermoTypeName); 

        // Print error message if package not found in the table
        //if (cstrIter == tablePtr->end())
        if(!ctorIter.good())
        {
            FatalErrorIn(Thermo::typeName + "::New")
                << "Unknown " << Thermo::typeName << " type " << nl
                << "thermoType" << thermoTypeDict << nl << nl
                << "Valid " << Thermo::typeName << " types are:" << nl << nl;

            // Get the list of all the suitable thermo packages available
            wordList validThermoTypeNames
            (
                tablePtr->sortedToc()
            );

            // Build a table of the thermo packages constituent parts
            // Note: row-0 contains the names of constituent parts
            List<wordList> validThermoTypeNameCmpts
            (
                validThermoTypeNames.size() + 1
            );

            validThermoTypeNameCmpts[0].setSize(nCmpt);
            forAll(validThermoTypeNameCmpts[0], j)
            {
                validThermoTypeNameCmpts[0][j] = cmptNames[j];
            }

            // Split the thermo package names into their constituent parts
            forAll(validThermoTypeNames, i)
            {
                validThermoTypeNameCmpts[i+1] =
                    Thermo::splitThermoName(validThermoTypeNames[i], nCmpt);
            }

            // Print the table of available packages
            // in terms of their constituent parts
            printTable(validThermoTypeNameCmpts, FatalError);

            FatalError<< exit(FatalError);
        }

        return ctorIter.val();
    }
    else
    {
        thermoTypeName = word(thermoDict.lookup("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

        //typename Table::iterator cstrIter = tablePtr->find(thermoTypeName);
        auto ctorIter = tablePtr->cfind(thermoTypeName);
        if (!ctorIter.good())
        {
            FatalErrorIn(Thermo::typeName + "::New")
                << "Unknown " << Thermo::typeName << " type "
                << thermoTypeName << nl << nl
                << "Valid " << Thermo::typeName << " types are:" << nl
                << tablePtr->sortedToc() << nl
                << exit(FatalError);
        }

        return ctorIter.val();
    }
}
*/


template<class Thermo, class ThermoConstructTable>
typename ThermoConstructTable::mapped_type
Foam::basic2Thermo::getThermoOrDie
(
    const dictionary& thermoTypeDict,
    ThermoConstructTable& thermoTable,
    const word& thermoTypeName,
    const wordList& cmptNames
)
{
    // Lookup the thermo package

    auto ctorIter = thermoTable.cfind(thermoTypeName);

    // Print error message if package not found in the table
    if (!ctorIter.good())
    {
        FatalIOErrorInLookup
        (
            thermoTypeDict,
            Thermo::typeName,
            word::null, // Suppress long name? Just output dictionary (above)
            thermoTable
        );

        basic2Thermo::printThermoNames
        (
            FatalIOError,
            cmptNames,
            thermoTable.sortedToc()
        ) << exit(FatalIOError);

        // return nullptr;
    }

    return ctorIter.val();
}


template<class Thermo, class ThermoConstructTable>
typename ThermoConstructTable::mapped_type
Foam::basic2Thermo::getThermoOrDie
(
    const dictionary& thermoDict,
    ThermoConstructTable& thermoTable
)
{
    const dictionary* dictptr = thermoDict.findDict("thermoType");

    if (dictptr)
    {
        const auto& thermoTypeDict = *dictptr;

        const wordList* cmptHeaderPtr = &(wordList::null());

        // Thermo package name, constructed from components
        const word thermoTypeName
        (
            basic2Thermo::makeThermoName(thermoTypeDict, cmptHeaderPtr)
        );

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        return getThermoOrDie<Thermo, ThermoConstructTable>
        (
            thermoTypeDict,
            thermoTable,
            thermoTypeName,
            *cmptHeaderPtr
        );
    }
    else
    {
        const word thermoTypeName(thermoDict.get<word>("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeName << endl;

        auto ctorIter = thermoTable.cfind(thermoTypeName);

        if (!ctorIter.good())
        {
            FatalIOErrorInLookup
            (
                thermoDict,
                Thermo::typeName,
                thermoTypeName,
                thermoTable
            ) << exit(FatalIOError);
        }

        return ctorIter.val();
    }
}

/*
template<class Thermo>
Foam::autoPtr<Thermo> Foam::basic2Thermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    //typename Thermo::fvMeshConstructorTable::iterator cstrIter =
    //    lookupThermo<Thermo, typename Thermo::fvMeshConstructorTable>
    //    (
    //        thermoDict,
    //        Thermo::fvMeshConstructorTablePtr_
    //    );

    //return autoPtr<Thermo>(cstrIter()(mesh, phaseName));

    auto* ctorPtr = lookupThermo<Thermo, typename Thermo::fvMeshConstructorTableType>
    (
        thermoDict,
        ((Thermo::fvMeshConstructorTablePtr_))
    );

    return autoPtr<Thermo>(ctorPtr(mesh, phaseName));

}
*/
/*
template<class Thermo>
Foam::autoPtr<Thermo> Foam::basic2Thermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    //typename Thermo::dictionaryConstructorTable::iterator cstrIter =
    //    lookupThermo<Thermo, typename Thermo::dictionaryConstructorTable>
    //    (
    //        dict,
    //        Thermo::dictionaryConstructorTablePtr_
    //    );
//
    //return autoPtr<Thermo>(cstrIter()(mesh, dict, phaseName));
    
    auto* ctorPtr = lookupThermo<Thermo, typename Thermo::dictionaryConstructorTableType>
    (
        dict,
        ((Thermo::dictionaryConstructorTablePtr_))
    );

    return autoPtr<Thermo>(ctorPtr(mesh, dict, phaseName));

}
    */
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::autoPtr<Thermo> Foam::basic2Thermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            phasePropertyName(dictName, phaseName),
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    auto* ctorPtr = getThermoOrDie<Thermo>
    (
        thermoDict,
        *(Thermo::fvMeshConstructorTablePtr_)
    );

    return autoPtr<Thermo>(ctorPtr(mesh, phaseName));
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basic2Thermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    auto* ctorPtr = getThermoOrDie<Thermo>
    (
        dict,
        *(Thermo::dictionaryConstructorTablePtr_)
    );

    return autoPtr<Thermo>(ctorPtr(mesh, dict, phaseName));
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basic2Thermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    IOdictionary thermoDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    auto* ctorPtr = getThermoOrDie<Thermo>
    (
        thermoDict,
        *(Thermo::fvMeshDictPhaseConstructorTablePtr_)
    );

    return autoPtr<Thermo>(ctorPtr(mesh, phaseName, dictName));
}



// ************************************************************************* //
