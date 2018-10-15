
# ifndef __MADELEINE_COMMONS_HPP__
# define __MADELEINE_COMMONS_HPP__

// ========================================================================== //
// INCLUDE                                                                    //
// ========================================================================== //

// Standard Template Library
# include <array>
# include <vector>

// bitpit
# include <bitpit_containers.hpp>
# include <bitpit_patchkernel.hpp>

// ========================================================================== //
// TYPES DEFINITIONS                                                          //
// ========================================================================== //

// C++ v11.0 arrays
typedef std::array<double, 1>       darray1;
typedef std::array<double, 2>       darray2;
typedef std::array<double, 3>       darray3;
typedef std::array<double, 5>       darray5;
typedef std::array<double, 7>       darray7;

typedef std::array<darray3,1>       darray31;
typedef std::array<darray3,2>       darray32;
typedef std::array<darray3,3>       darray33;
typedef std::array<darray3,5>       darray35;
typedef std::array<darray3,7>       darray37;
typedef std::array<darray5,5>       darray55;

# endif
