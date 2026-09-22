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

#include "Cloud.H"
#include "Time.H"
#include "IOPosition.H"

#include <atomic>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParticleType>
Foam::word Foam::Cloud<ParticleType>::cloudPropertiesName("cloudProperties");


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class ParticleType>
void Foam::Cloud<ParticleType>::readCloudUniformProperties()
{
    IOobject dictObj
    (
        cloudPropertiesName,
        time().timeName(),
        "uniform"/cloud::prefix/name(),
        db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );

    //if (dictObj.headerOk())
    if (dictObj.typeHeaderOk<IOList<label>>(false)) // NEW VINCENT
    {
        const IOdictionary uniformPropsDict(dictObj);

        const word procName("processor" + Foam::name(Pstream::myProcNo()));
        if (uniformPropsDict.found(procName))
        {
            uniformPropsDict.subDict(procName).lookup("particleCount")
                >> ParticleType::particleCount_;
        }
    }
    else
    {
        ParticleType::particleCount_ = 0;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeCloudUniformProperties() const
{
    IOdictionary uniformPropsDict
    (
        IOobject
        (
            cloudPropertiesName,
            time().timeName(),
            "uniform"/cloud::prefix/name(),
            db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelList np(Pstream::nProcs(), 0);
    np[Pstream::myProcNo()] = ParticleType::particleCount_;

    Pstream::listCombineGather(np, maxEqOp<label>());
    Pstream::listCombineScatter(np);

    forAll(np, i)
    {
        word procName("processor" + Foam::name(i));
        uniformPropsDict.add(procName, dictionary());
        uniformPropsDict.subDict(procName).add("particleCount", np[i]);
    }

    uniformPropsDict.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        time().writeCompression()
    );
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::initCloud(const bool checkClass)
{
    readCloudUniformProperties();

    IOPosition<Cloud<ParticleType> > ioP(*this);

    if (ioP.headerOk())
    {
        ioP.readData(*this, checkClass);
        ioP.close();

        if (this->size())
        {
            readFields();
        }
    }
    else
    {
        if (debug)
        {
            Pout<< "Cannot read particle positions file:" << nl
                << "    " << ioP.objectPath() << nl
                << "Assuming the initial cloud contains 0 particles." << endl;
        }
    }

    initCloudPostRead();
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::initCloudPostRead(const bool reLocate)
{
    // Ask for the tetBasePtIs to trigger all processors to build
    // them, otherwise, if some processors have no particles then
    // there is a comms mismatch.
    polyMesh_.tetBasePtIs();

    if (!reLocate)
    {
        // Tet data already computed per particle (parallel filtered read)
        return;
    }

    scalar lostParticles = 0; // NEW VINCENT
    scalar totParticles = 0; // NEW VINCENT

    forAllIter(typename Cloud<ParticleType>, *this, pIter)
    {
        ParticleType& p = pIter();
        ++totParticles;

        //p.initCellFacePt();
        // Mass will not be conserved but might be OK if sampling is not
        // yet enabled (because refinement steps are still in progress)
        if (p.initCellFacePtOrDeleteLostParticle()) // NEW VINCENT
        {
            ++lostParticles;
            deleteParticle(p);
        }
    }

    if (lostParticles > 0) // NEW VINCENT
    {
        Info<< "Lost particles deleted due to change in topology:" << tab
            << lostParticles << "/" << totParticles << endl;
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::readCloudFiltered
(
    const labelUList& cellOwner,
    const label filterRank
)
{
    initialReadKeep_.clear();
    initialReadNFull_ = 0;

    readCloudUniformProperties();

    IOPosition<Cloud<ParticleType> > ioP(*this);

    if (ioP.headerOk())
    {
        const fileName posFile = ioP.objectPath();
        off_t fileSize = 0;
        {
            const int fd = ::open(posFile.c_str(), O_RDONLY);
            if (fd >= 0)
            {
                struct stat st;
                if (::fstat(fd, &st) == 0)
                {
                    fileSize = st.st_size;
                }
                ::close(fd);
            }
        }

        label nThreads = 1;
        #ifdef _OPENMP
        nThreads = omp_get_max_threads();
        #endif

        const bool useParallel =
        (
            nThreads > 1
         && fileSize > (off_t(100)<<20)   // small files: serial read is faster
        );

        bool parallelOk = false;
        if (useParallel)
        {
            parallelOk =
                parallelFilteredRead
                (
                    ioP,
                    posFile,
                    fileSize,
                    cellOwner,
                    filterRank,
                    initialReadKeep_,
                    initialReadNFull_
                );
        }

        if (!parallelOk)
        {
            // checkClass=false: same as the dsmcCloud constructor path
            // (Cloud initCloud(false)); the positions header class is
            // "Cloud" while IOPosition's inherited typeName would compare
            // as "regIOobject".
            ioP.readDataFiltered
            (
                *this,
                false,
                cellOwner,
                filterRank,
                initialReadKeep_,
                initialReadNFull_
            );
            ioP.close();
            initCloudPostRead(true);
        }
        else
        {
            // tet data already computed per particle in the parallel read
            initCloudPostRead(false);
        }
    }
    else
    {
        if (debug)
        {
            Pout<< "Cannot read particle positions file:" << nl
                << "    " << ioP.objectPath() << nl
                << "Assuming the initial cloud contains 0 particles." << endl;
        }
        initCloudPostRead(false);
    }
}


template<class ParticleType>
bool Foam::Cloud<ParticleType>::parallelFilteredRead
(
    IOPosition<Cloud<ParticleType>>& ioP,
    const fileName& posFile,
    const off_t fileSize,
    const labelUList& cellOwner,
    const label filterRank,
    labelList& keep,
    label& nFull
)
{
    // ---- header, entry count, entries region start ----
    off_t entriesStart = 0;
    nFull = ioP.readHeader(entriesStart);
    ioP.close();

    if (nFull < 1)
    {
        keep.clear();
        return true;    // nothing to read
    }

    // ---- mmap the positions file (page cache shared across ranks) ----
    const int fd = ::open(posFile.c_str(), O_RDONLY);
    if (fd < 0)
    {
        FatalIOErrorInFunction(posFile)
            << "cannot open " << posFile << exit(FatalIOError);
    }
    const char* base = static_cast<const char*>
    (
        ::mmap(nullptr, size_t(fileSize), PROT_READ, MAP_PRIVATE, fd, 0)
    );
    ::close(fd);
    if (base == MAP_FAILED)
    {
        FatalIOErrorInFunction(posFile)
            << "cannot mmap " << posFile << exit(FatalIOError);
    }

    // ---- chunk plan: nThreads chunks of ceil(nFull/nThreads) entries ----
    label nThreads = 1;
    #ifdef _OPENMP
    nThreads = omp_get_max_threads();
    #endif
    nThreads = max(1, min(nThreads, nFull));
    const label chunkEntries = (nFull + nThreads - 1)/nThreads;

    List<off_t> byteStart(nThreads + 1, off_t(-1));
    List<label> idxStart(nThreads + 1, label(-1));
    byteStart[0] = entriesStart;
    idxStart[0] = 0;
    for (label t = 1; t <= nThreads; ++t)
    {
        idxStart[t] = min(t*chunkEntries, nFull);
    }

    // ---- single-pass scan: byte offset of each chunk's first entry ----
    {
        const char* p = base + entriesStart;
        const char* fileEnd = base + fileSize;
        label entryI = 0;
        label nextTarget = 1;
        bool inEntries = true;

        while (p < fileEnd && nextTarget <= nThreads)
        {
            const char* nl =
                static_cast<const char*>(memchr(p, '\n', fileEnd - p));
            const char* lineEnd = (nl ? nl : fileEnd);

            if (*p == '(')
            {
                if (entryI == idxStart[nextTarget])
                {
                    byteStart[nextTarget] = p - base;
                    ++nextTarget;
                    if (nextTarget > nThreads)
                    {
                        break;
                    }
                }
                ++entryI;
            }
            else if (*p == ')')
            {
                inEntries = false;
                break;
            }

            p = lineEnd + 1;
        }

        if (inEntries || entryI != nFull)
        {
            ::munmap(const_cast<char*>(base), size_t(fileSize));
            FatalIOErrorInFunction(posFile)
                << "positions file entry scan found " << entryI
                << " entries, expected " << nFull << exit(FatalIOError);
        }
        byteStart[nThreads] = (p - base) + 1;
        idxStart[nThreads] = nFull;
    }

    // ---- warm-up mesh location caches (single-threaded) ----
    polyMesh_.tetBasePtIs();

    Info<< "Cloud filtered read: parsing " << nFull << " entries with "
        << nThreads << " threads" << endl;

    // ---- parallel parse + filter ----
    List<DynamicList<ParticleType*>> keptLocal(nThreads);
    List<DynamicList<label>> keepLocal(nThreads);
    label parsedTotal = 0;

    #ifdef _OPENMP
    #pragma omp parallel num_threads(nThreads)
    #endif
    {
        label t = 0;
        #ifdef _OPENMP
        t = omp_get_thread_num();
        #endif
        DynamicList<ParticleType*>& kept = keptLocal[t];
        DynamicList<label>& kp = keepLocal[t];

        const char* p = base + byteStart[t];
        const char* end = base + byteStart[t + 1];
        label gi = idxStart[t];

        while (p < end)
        {
            const char* nl =
                static_cast<const char*>(memchr(p, '\n', end - p));
            const char* lineEnd = (nl ? nl : end);

            if (*p == '(')
            {
                // entry format written by IOPosition: "(x y z) cellI"
                char* q = const_cast<char*>(p + 1);
                const double x = ::strtod(q, &q);
                const double y = ::strtod(q, &q);
                const double z = ::strtod(q, &q);
                while (q < lineEnd && (*q == ' ' || *q == ')'))
                {
                    ++q;
                }
                const label cellI = label(::strtol(q, &q, 10));

                if (cellI < 0 || cellI >= cellOwner.size())
                {
                    FatalErrorInFunction
                        << "Entry " << gi << " has out-of-range cell index "
                        << cellI << exit(FatalError);
                }

                if (cellOwner[cellI] == filterRank)
                {
                    ParticleType* pPtr =
                        new ParticleType
                        (
                            polyMesh_,
                            vector(x, y, z),
                            cellI,
                            true    // locate via initCellFacePtOrDeleteLostParticle
                        );

                    if (pPtr->initCellFacePtOrDeleteLostParticle())
                    {
                        delete pPtr;    // lost particle
                    }
                    else
                    {
                        kept.append(pPtr);
                        kp.append(gi);
                    }
                }
                // else: non-owned, skipped without construction cost

                ++gi;

                label done = 0;
                #ifdef _OPENMP
                #pragma omp atomic capture
                #endif
                done = ++parsedTotal;

                if ((done % 20000000) == 0)
                {
                    #pragma omp critical(progressPrint)
                    {
                        Info<< "Cloud filtered read: parsed " << done
                            << " / " << nFull << endl;
                    }
                }
            }
            else if (*p == ')')
            {
                break;
            }

            p = lineEnd + 1;
        }
    }

    ::munmap(const_cast<char*>(base), size_t(fileSize));

    // ---- serial merge in chunk order (deterministic cloud layout) ----
    for (label t = 0; t < nThreads; ++t)
    {
        forAll(keptLocal[t], i)
        {
            this->append(keptLocal[t][i]);
            keep.append(keepLocal[t][i]);
        }
        keptLocal[t].clear();
        keepLocal[t].clear();
    }

    Info<< "Cloud filtered read: kept " << this->size()
        << " / " << nFull << " particles" << endl;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const bool checkClass
)
:
    cloud(pMesh),
    polyMesh_(pMesh),
    labels_(),
    nTrackingRescues_(),
    cellWallFacesPtr_(),
    openmpMoveMeshDataReady_(false),
    initialReadKeep_(),
    initialReadNFull_(0)
{
    checkPatches();

    initCloud(checkClass);
}


template<class ParticleType>
Foam::Cloud<ParticleType>::Cloud
(
    const polyMesh& pMesh,
    const word& cloudName,
    const bool checkClass,
    const bool readPositions
)
:
    cloud(pMesh, cloudName),
    polyMesh_(pMesh),
    labels_(),
    nTrackingRescues_(),
    cellWallFacesPtr_(),
    openmpMoveMeshDataReady_(false),
    initialReadKeep_(),
    initialReadNFull_(0)
{
    checkPatches();

    if (readPositions)
    {
        initCloud(checkClass);
    }
    else
    {
        // Deferred read: uniform properties only; positions and fields
        // are read later via readCloudFiltered + readFieldsFiltered.
        readCloudUniformProperties();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParticleType>
Foam::IOobject Foam::Cloud<ParticleType>::fieldIOobject
(
    const word& fieldName,
    const IOobject::readOption r
) const
{
    return IOobject
    (
        fieldName,
        time().timeName(),
        *this,
        r,
        IOobject::NO_WRITE,
        false
    );
}


template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldIOobject
(
    const Cloud<ParticleType>& c,
    const IOField<DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorIn
        (
            "void Cloud<ParticleType>::checkFieldIOobject"
            "(const Cloud<ParticleType>&, const IOField<DataType>&) const"
        )   << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
template<class DataType>
void Foam::Cloud<ParticleType>::checkFieldFieldIOobject
(
    const Cloud<ParticleType>& c,
    const CompactIOField<Field<DataType>, DataType>& data
) const
{
    if (data.size() != c.size())
    {
        FatalErrorIn
        (
            "void Cloud<ParticleType>::checkFieldFieldIOobject"
            "("
                "const Cloud<ParticleType>&, "
                "const CompactIOField<Field<DataType>, DataType>&"
            ") const"
        )   << "Size of " << data.name()
            << " field " << data.size()
            << " does not match the number of particles " << c.size()
            << abort(FatalError);
    }
}


template<class ParticleType>
void Foam::Cloud<ParticleType>::readFields()
{}


template<class ParticleType>
void Foam::Cloud<ParticleType>::writeFields() const
{
    if (this->size())
    {
        ParticleType::writeFields(*this);
    }
}


template<class ParticleType>
bool Foam::Cloud<ParticleType>::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    writeCloudUniformProperties();

    if (this->size())
    {
        writeFields();
        return cloud::writeObject(fmt, ver, cmp);
    }
    else
    {
        return true;
    }
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

template<class ParticleType>
Foam::Ostream& Foam::operator<<(Ostream& os, const Cloud<ParticleType>& pc)
{
    pc.writeData(os);

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const Cloud<ParticleType>&)");

    return os;
}


// ************************************************************************* //
