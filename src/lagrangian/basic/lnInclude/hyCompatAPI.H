#ifndef hyCompatAPI_H
#define hyCompatAPI_H

#include "dictionary.H"
#include "label.H"
#include "scalar.H"
#include <type_traits>

namespace Foam
{
namespace hyCompat
{

class lookupResult
{
    ITstream& is_;

public:

    explicit lookupResult(ITstream& is)
    :
        is_(is)
    {}

    operator ITstream&() const
    {
        return const_cast<ITstream&>(is_);
    }

    template<class T>
    operator T() const
    {
        if constexpr (std::is_arithmetic_v<T>)
        {
            T value{};
            is_ >> value;
            return value;
        }
        else
        {
            return T(is_);
        }
    }
};


inline lookupResult lookup
(
    const dictionary& dict,
    const word& keyword,
    enum keyType::option matchOpt = keyType::REGEX
)
{
    return lookupResult(dict.lookupEntry(keyword, matchOpt).stream());
}



template<class T>
inline const T& lookup(const T& value)
{
    return value;
}
inline const entry& lookupEntry
(
    const dictionary& dict,
    const word& keyword,
    enum keyType::option matchOpt = keyType::REGEX
)
{
    return dict.lookupEntry(keyword, matchOpt);
}


template<class ListType, class T>
inline label indexOf(const ListType& list, const T& value)
{
    return list.find(value);
}


inline scalar toScalar(const lookupResult& is)
{
    ITstream& stream = is;
    scalar value = Zero;
    stream >> value;
    return value;
}


template<class StreamType, class = std::enable_if_t<!std::is_same_v<std::decay_t<StreamType>, lookupResult>>>
inline scalar toScalar(StreamType&& is)
{
    scalar value = Zero;
    is >> value;
    return value;
}


inline scalar toScalar(const scalar value)
{
    return value;
}


inline scalar toScalar
(
    const dictionary& dict,
    const word& keyword,
    enum keyType::option matchOpt = keyType::REGEX
)
{
    return dict.get<scalar>(keyword, matchOpt);
}


inline label toLabel(const lookupResult& is)
{
    ITstream& stream = is;
    label value = 0;
    stream >> value;
    return value;
}


template<class StreamType, class = std::enable_if_t<!std::is_same_v<std::decay_t<StreamType>, lookupResult>>>
inline label toLabel(StreamType&& is)
{
    label value = 0;
    is >> value;
    return value;
}


inline label toLabel(const label value)
{
    return value;
}


inline label toLabel
(
    const dictionary& dict,
    const word& keyword,
    enum keyType::option matchOpt = keyType::REGEX
)
{
    return dict.get<label>(keyword, matchOpt);
}

} // End namespace hyCompat
} // End namespace Foam

#endif

