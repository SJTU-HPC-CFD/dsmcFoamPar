#include <metis.h>
#include <mpi.h>

#include <cfenv>
#include <csignal>
#include <cstdlib>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct Args
{
    std::string ownerPath;
    std::string neighbourPath;
    int cells = -1;
    int parts = -1;
    bool enableFpe = false;
    bool verbose = false;
};


void signalHandler(int sig)
{
    std::cerr << "[metis-minimal] caught signal " << sig << std::endl;
    std::_Exit(128 + sig);
}


Args parseArgs(int argc, char** argv)
{
    Args args;

    for (int i = 1; i < argc; ++i)
    {
        const std::string token(argv[i]);

        if (token == "--owner" && i + 1 < argc)
        {
            args.ownerPath = argv[++i];
        }
        else if (token == "--neighbour" && i + 1 < argc)
        {
            args.neighbourPath = argv[++i];
        }
        else if (token == "--cells" && i + 1 < argc)
        {
            args.cells = std::atoi(argv[++i]);
        }
        else if (token == "--parts" && i + 1 < argc)
        {
            args.parts = std::atoi(argv[++i]);
        }
        else if (token == "--fpe")
        {
            args.enableFpe = true;
        }
        else if (token == "--verbose")
        {
            args.verbose = true;
        }
        else
        {
            throw std::runtime_error("unknown or incomplete option: " + token);
        }
    }

    if (args.ownerPath.empty() || args.neighbourPath.empty())
    {
        throw std::runtime_error("both --owner and --neighbour are required");
    }

    if (args.parts <= 0)
    {
        throw std::runtime_error("--parts must be > 0");
    }

    return args;
}


std::vector<int> readLabelList(const std::string& path)
{
    std::ifstream in(path.c_str());

    if (!in)
    {
        throw std::runtime_error("cannot open " + path);
    }

    std::string line;

    for (int i = 0; i < 10; ++i)
    {
        if (!std::getline(in, line))
        {
            throw std::runtime_error("unexpected EOF while reading header: " + path);
        }
    }

    int expected = 0;
    {
        std::istringstream iss(line);
        iss >> expected;
    }

    if (expected <= 0)
    {
        throw std::runtime_error("invalid list size in " + path + ": " + line);
    }

    if (!std::getline(in, line))
    {
        throw std::runtime_error("missing opening parenthesis in " + path);
    }

    std::vector<int> values;
    values.reserve(expected);

    while (static_cast<int>(values.size()) < expected && std::getline(in, line))
    {
        if (line.find(')') != std::string::npos)
        {
            break;
        }

        std::istringstream iss(line);
        int value = 0;
        while (iss >> value)
        {
            values.push_back(value);
        }
    }

    if (static_cast<int>(values.size()) != expected)
    {
        std::ostringstream oss;
        oss << "size mismatch in " << path
            << ": expected " << expected
            << ", got " << values.size();
        throw std::runtime_error(oss.str());
    }

    return values;
}


struct CellStats
{
    int minOwner = std::numeric_limits<int>::max();
    int maxOwner = std::numeric_limits<int>::min();
    int minNeighbour = std::numeric_limits<int>::max();
    int maxNeighbour = std::numeric_limits<int>::min();
};


CellStats computeCellStats
(
    const std::vector<int>& owner,
    const std::vector<int>& neighbour
)
{
    CellStats stats;

    for (const int value : owner)
    {
        if (value < stats.minOwner) stats.minOwner = value;
        if (value > stats.maxOwner) stats.maxOwner = value;
    }

    for (const int value : neighbour)
    {
        if (value < stats.minNeighbour) stats.minNeighbour = value;
        if (value > stats.maxNeighbour) stats.maxNeighbour = value;
    }

    return stats;
}


void enableFpeIfRequested(bool enable)
{
    if (!enable)
    {
        return;
    }

    std::signal(SIGFPE, signalHandler);

    #if defined(__linux__)
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    #endif
}

} // namespace


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int exitCode = 0;

    try
    {
        const Args args = parseArgs(argc, argv);
        enableFpeIfRequested(args.enableFpe);

        if (rank == 0)
        {
            std::cout
                << "[metis-minimal] rank0 start"
                << " owner=" << args.ownerPath
                << " neighbour=" << args.neighbourPath
                << " cells=" << args.cells
                << " parts=" << args.parts
                << " mpiSize=" << size
                << " fpe=" << (args.enableFpe ? 1 : 0)
                << std::endl;

            const std::vector<int> owner = readLabelList(args.ownerPath);
            const std::vector<int> neighbour = readLabelList(args.neighbourPath);

            if (owner.size() < neighbour.size())
            {
                throw std::runtime_error("owner size < neighbour size");
            }

            const CellStats stats = computeCellStats(owner, neighbour);
            const int inferredCells =
                std::max(stats.maxOwner, stats.maxNeighbour) + 1;
            const int nCells =
                (args.cells > 0 ? args.cells : inferredCells);
            const int nIntFaces = static_cast<int>(neighbour.size());

            std::cout
                << "[metis-minimal] inferred cells=" << inferredCells
                << " minOwner=" << stats.minOwner
                << " maxOwner=" << stats.maxOwner
                << " minNeighbour=" << stats.minNeighbour
                << " maxNeighbour=" << stats.maxNeighbour;

            if (args.cells > 0)
            {
                std::cout << " overrideCells=" << args.cells;
            }

            std::cout << std::endl;

            std::vector<int> degree(nCells, 0);
            for (int faceI = 0; faceI < nIntFaces; ++faceI)
            {
                const int own = owner[faceI];
                const int nei = neighbour[faceI];

                if (own < 0 || own >= nCells || nei < 0 || nei >= nCells)
                {
                    std::ostringstream oss;
                    oss << "face " << faceI
                        << " has out-of-range cells: own=" << own
                        << " nei=" << nei
                        << " nCells=" << nCells;
                    throw std::runtime_error(oss.str());
                }

                ++degree[own];
                ++degree[nei];
            }

            std::vector<idx_t> xadj(nCells + 1, 0);
            for (int cellI = 0; cellI < nCells; ++cellI)
            {
                xadj[cellI + 1] = xadj[cellI] + degree[cellI];
            }

            std::vector<idx_t> adjncy(static_cast<std::size_t>(xadj[nCells]));
            std::vector<int> offset(nCells, 0);
            for (int cellI = 0; cellI < nCells; ++cellI)
            {
                offset[cellI] = static_cast<int>(xadj[cellI]);
            }

            for (int faceI = 0; faceI < nIntFaces; ++faceI)
            {
                const int own = owner[faceI];
                const int nei = neighbour[faceI];
                adjncy[offset[own]++] = nei;
                adjncy[offset[nei]++] = own;
            }

            idx_t nvtxs = nCells;
            idx_t ncon = 1;
            idx_t nparts = args.parts;
            idx_t edgecut = 0;
            std::vector<idx_t> part(nCells, 0);
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions(options);
            options[METIS_OPTION_SEED] = 42;

            const int rc = METIS_PartGraphKway
            (
                &nvtxs,
                &ncon,
                xadj.data(),
                adjncy.data(),
                nullptr,
                nullptr,
                nullptr,
                &nparts,
                nullptr,
                nullptr,
                options,
                &edgecut,
                part.data()
            );

            long long checksum = 0;
            idx_t minPart = part.empty() ? idx_t(-1) : part[0];
            idx_t maxPart = part.empty() ? idx_t(-1) : part[0];
            for (std::size_t i = 0; i < part.size(); ++i)
            {
                checksum += static_cast<long long>(part[i]) * static_cast<long long>(i + 1);
                if (part[i] < minPart) minPart = part[i];
                if (part[i] > maxPart) maxPart = part[i];
            }

            std::cout
                << "[metis-minimal] METIS_PartGraphKway rc=" << rc
                << " edgecut=" << edgecut
                << " minPart=" << minPart
                << " maxPart=" << maxPart
                << " checksum=" << checksum
                << std::endl;

            if (args.verbose)
            {
                std::cout << "[metis-minimal] first 16 parts:";
                const std::size_t limit = std::min<std::size_t>(16, part.size());
                for (std::size_t i = 0; i < limit; ++i)
                {
                    std::cout << ' ' << part[i];
                }
                std::cout << std::endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }
    catch (const std::exception& e)
    {
        std::cerr << "[metis-minimal] rank " << rank
                  << " exception: " << e.what() << std::endl;
        exitCode = 2;
    }
    catch (...)
    {
        std::cerr << "[metis-minimal] rank " << rank
                  << " unknown exception" << std::endl;
        exitCode = 3;
    }

    int globalExit = 0;
    MPI_Allreduce(&exitCode, &globalExit, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Finalize();
    return globalExit;
}
