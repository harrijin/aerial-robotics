#include "typedefs.h"

#ifndef __DEFINITIONS_H
#define __DEFINITIONS_H

//======= General
constexpr s32 SEC_PER_MINUTE = 60;
constexpr s32 MINUTE_PER_HOUR = 60;
constexpr s32 SEC_PER_HOUR  = 3600;
constexpr s32 SEC_PER_DAY   = 86400;
constexpr s32 SEC_PER_WEEK  = 604800;
constexpr s32 DAYS_PER_WEEK = 7;
constexpr s32 HOURS_PER_DAY = 24;

//======= Enumeration types meant to be widely accessible
enum ReturnValue {
  RETVAL_SUCCESS = 0,
  RETVAL_FAILURE = -1
};

//======= Immutable constants
constexpr f64 PI = 3.14159265358979;
constexpr f64 TWOPI = 2*PI;
constexpr s32 SPEED_LIGHT_MPS = 299792458;
constexpr f64 DEG_TO_RAD = PI/180.0;
constexpr f64 RAD_TO_DEG = 180.0/PI; 

//======= WGS84 constants
extern const f64 sqrt_muearth_WGS84;
constexpr f64 AA_WGS84 = 6378137.00000;		   // meters
constexpr f64 BB_WGS84 = 6356752.31425;		   // meters
constexpr f64 OmegaE_WGS84 = 7.2921151467e-5;	   // radians/second
constexpr f64 e_WGS84 = 0.0818191908334158;        // eccentricity
constexpr f64 esquare_WGS84 = 0.00669437998863492; // eccenctricity squared

#endif

