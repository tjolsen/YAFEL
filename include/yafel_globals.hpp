#ifndef YAFEL_GLOBALS_HPP
#define YAFEL_GLOBALS_HPP

#include<cstddef>

/**
 * \file
 *
 * This header defines any global macros, types, etc that are used throughout the library
*/

// project namespace open/close macros
#define YAFEL_NAMESPACE_OPEN namespace yafel {
#define YAFEL_NAMESPACE_CLOSE }

namespace yafel {

using size_type = int;


#ifndef YAFEL_ALWAYS_INLINE
#define YAFEL_ALWAYS_INLINE __attribute__((always_inline)) inline
#endif

#ifndef YAFEL_NEVER_INLINE
#define YAFEL_NEVER_INLINE __attribute__((noinline))
#endif


} // end namespace yafel

#endif // end #ifndef YAFEL_GLOBALS_HPP
