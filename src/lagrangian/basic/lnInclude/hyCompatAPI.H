#ifndef hyCompatAPI_H
#define hyCompatAPI_H

#include "dictionary.H"
#include "label.H"
#include "scalar.H"
#include "word.H"
#include "UList.H"
#include "stringOps.H"
#include "token.H"
#include "ITstream.H"

namespace Foam
{
namespace hyCompat
{

template<class DictType>
inline auto lookup(const DictType& dict, const word& key)
-> decltype(dict.lookup(key))
{
    return dict.lookup(key);
}

template<class T>
inline T& lookup(T& value)
{
    return value;
}

template<class T>
inline const T& lookup(const T& value)
{
    return value;
}

template<class DictType>
inline scalar toScalar(const DictType& dict, const word& key)
{
    return dict.template get<scalar>(key);
}

inline scalar toScalar(Istream& is)
{
    return readScalar(is);
}

inline scalar toScalar(const token& tok)
{
    return tok.number();
}

template<class DictType>
inline label toLabel(const DictType& dict, const word& key)
{
    return dict.template get<label>(key);
}

inline label toLabel(Istream& is)
{
    return readLabel(is);
}

inline label toLabel(const token& tok)
{
    return tok.labelToken();
}

template<class ListType, class T>
inline label indexOf(const ListType& list, const T& value)
{
    return findIndex(list, value);
}

}
}

#endif
