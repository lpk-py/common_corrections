.\" @(#)manual.l;
.\"===========================================================================*
.\" Program for xxxx    |                raydyntrace                 |      manual page |
.\"===========================================================================*
.\"                *=======================================*
.\"================|     Manual page header definition     |===================
.\"                *=======================================*
.TH PROGRAM 1
.\"                *=======================================*
.\"================|     Command name and description      |===================
.\"                *=======================================*
.SH NAME
raydyntrace \- dynamic ray tracing in a spherical Earth
.\"                *=======================================*
.\"================|           Command synopsis            |===================
.\"                *=======================================*
.\" a short synopsis of the command syntax, including its options and arguments
.SH SYNOPSIS
.B raydyntrace
.br
.nf
.\"================|         Command description           |===================
.\"                *=======================================*
.SH DESCRIPTION
.B raydyntrace
does dynamic raytracing in a spherically symmetric Earth. It computes
at the same time topographic and ellipticity corrections for arbitrary ray definitions.
Some input is asked from the screen. It can be used in two ways: to compute rays
among individual source-station paths, in which case the output file is compatible
with the output of program raydata.f, used for the computation of finite-frequency
kernels. Or it can generate a table for more general use.

.B
Using raydyntrace to compute a ray table
.br

Calvet and Chevrot (PEPI 153:21-31, 2005) compute a table of rays, with travel times
and geometrical spreading factors, to avoid the paraxial approximation for core phases 
arriving near 180 degrees epicentral distance (when the arrival is sensitive to the 
structure over all azimuths and not just in the neighbourhood of the geometrical ray 
path). The same may be neede for phases like PP or SS near the antipode. For this
purpose, Raydyntrace can also be run in "table generation mode" by setting the tableflag > 0
- the last of the screne input parameters - larger than 0. The program will then attempt 
to create a table with rays roughly spaced at intervals of 1/tableflag degrees (e.g.
if you set tableflag equal to 4, the rays will be spaced approximately 0.25 degree
in epicentral distance). You should define the source depth correctly in the ray definition
file (in contrast to the regular use of raydyntrace when the source depth is adapted to
that listed in the raypaths file).
For a discussion of the issue of ray sensitivity near the antipode see
also Tian et al., J. Comp. Phys., 226:2271-2288, 2007.


.\"                *=======================================*
.\"================|         Files used on input           |===================
.\"                *=======================================*
.\" names of files that raydyntrace creates
.SH INPUT FILES
Before you run the program, you need:
.br
a. A model file (e.g. IASP91.txt)
.br
b. A ray definition file (e.g. Pdef)
.br
c. A file containing source- and receiver information for each raypath (e.g. Pshallow)
unless you use the tableflag>0 option to create a table
.br
d. Crustal model CRUST2.0 (see mahi.ucsd.edu/Gabi/rem.dir/crust/crust2.html)
unless you use the tableflag>0 option to create a table
.br

These files have free formats:

.I
a. Model file

       Line 1: ignored
.br
       Line 2: ignored
.br
       Line 3: n (# of nodes), nic (top node in inner core), noc (top
.br
               node in outer core), nmoh (bottom[!] node of the crust)
.br
       Lines 4ff: radius, density, Vp, Vs
.br

.br
The nodes are ordered with INCREASING radius.
Discontinuities are indicated by two subsequent nodes
with the same radius. Nodes closer than 0.01 km are collapsed
into a discontinuity. However, because we use spline interpolation
between nodes, take care with node spacing (we experienced
a ray caught in a spline-created low velocity zone when we replaced
an upper mantle discontinuity by a 10 km gradient - if you wish to 
do this, space your nodes such that the spline does what you intended).
Model units are assumed to be MKS, unless the first density
read is less than 100, in which case km and km/s units are
assumed. This is the same model format used in several normal
mode codes with the exception of the addition of the moho layer
in line 3.

.I
b. Ray definition file
.br

The ray is defined by segments and their starting radii.
The first (source) one specifies: radius, up/down, P/S.
For each subsequent segment: radius,refl/trans/turn,P/S
and for the last point of the ray: (radius,5,0)). The source
depth can be arbitrary (e.g. 6371) and will be adapted to
reflect the real source depth for each ray in the source/receiver
file.
Note that a turning point separates two different segments
(the radius is then interpreted as the MINIMUM radius for 
the turning point).

      Line 1: Descriptive line (max 72 characters)
.br
      Line 2: Ray mnemonic (e.g. ScP)
.br
      Line 3ff: Segment lines (see below)
.br

Segments are separated by either reflection, transmission, or turning
points. If there is no wave conversion at a layer, one does not need
to specify a transmission.
The last segment is not a real segment but indicates the point where
the previous segment should end (usually the surface of the Earth).
Each segment line contains: radius (km), ray direction, wave type:
.br
- the radius is the starting radius of the segment, or the minimum
  possible radius in case of a turning point,
.br
- the direction is coded as follows:1=down from source,0=up from source,
  2=turning point,3=refl, 4=transmission, 5=last point of ray
.br
- the wave type is either 1 (P) or 2 (S). For historical reasons
  the last wave type may be 0, in which it will be reset to the
  same type as the last segment.


.I
c. Source/receiver file (if using option tableflag=0)
.br

This file has the coordinates in degrees and the depth/receiver
elevation in km. For each ray there is a line:
Date, Source latitude, longitude, depth, Receiver latitude,
longitude, elevation, Receiver code, kpole.
Date is an integer up to i8, which the user can use to identify the event.
Latitudes and longitudes are in degrees. The receiver code or name is
a stationcode of up to 16 characters, although only 8 are preserved
in the output file. Finally, kpole gives the number of polar crossings
for the ray. E.g., for a shortest epicentral distance of 150 deg and kpole
equal to 0, rays like PPP will be found that travel 150 degrees. But if kpole
equals one, the effective distance is assumed to be 210 deg.

When running the program, the user is asked for the names of these three
files as well as for:
.br
- the maximum number of arrivals, if more then one, to include in the output,
.br
- the maximum delay that a later arrival may have to be included in the output.
.br
- the value of tableflag (see below)
.br


.\"                *=======================================*
.\"================|         Files used or created         |===================
.\"                *=======================================*
.\" names of files that raydyntrace creates
.SH FILES
.br

A summary of the computations is provided in the file out.raydyntrace.xxxx
if xxxx is the name of the file with raypaths.

The most important output file is the file with all ray information, including
the Hessian matrix (second derivatives of travel time when one moves away
from the minimum travel time path), summed for forward- and backwards travelling
ray. The name of this file is raydyntrace.xxxx if xxxx is the name of the file with
raypaths. 

.B
If you specified tableflag=0 (input from raypath file):
.br

The first few lines repeat the ray definition with comments. Directly following are,
for each ray:
.br
  Line 1ff: segment number, ray direction code, latitude, longitude, radius, crustal
          correction for CRUST2.0, elevation correction [one line per segment]
  (repeated for each segment)
.br
  Line 2: idate,source lat,lon,radius, receiver lat,lon,elevation,station code
.br
  Line 3: arrival number, number of ray nodes, 0(zero), travel time (s, uncorrected),
          ellipticity correction, t*, Q, slowness (rad/s)
.br
  Line 4ff: r,i,phi,seismic velocity at the ray node, 1/Q, H_11, H_22
           [one line per node]
.br
In 1,longitude and latitude of the segment nodes are in degrees, corrections are in
second and have to be added to the computed travel time for the 1D Earth to get
the travel time in an elliptical Earth with a variable crust.
.br
In 2, units are again degrees, km.
.br
In 3, the arrivals are numbered in order of increasing travel time. The zero is
for compatibility and can be ignored. The ellipticity correction and t* are in
seconds.
.br
In 4, r is in km, the ray angle i and epicentral distance phi are in radians,
the seismic velocity is in km/s, the Hessians H_11 and H_22 are in sec/km^2.
H_11 and H_22 are the diagonal elements of the sum of the forward- and backward 
Hessians M' and M" in Dahlen et al., 2000 (see also Nolet, 2008 or Tian et al., 2007).
.br

.B
If you specify tableflag>0 (generating a table for interpolation):

The format of the output file differs in case tableflag>0 (table generation mode).
The first few lines still repeat the ray definition, followed by, for every ray:
.br
  Line 1: a header line with pertinent information (ray nr, # of nodes, p, starting 
          angle,Rrs,Vs,Vr,Trs)
.br
  Line 2: a line with column headers
.br
  Line 3ff: r(km), i(rad), Delta(rad), T(sec), c(km/s), q0=1/Qs, h11 and h22 (sec/km^2),
          Rxs and Rxr (km) and length along ray (km).
.br
(NOTE: The version of Raydyntrace previous to Dec 2008 printed only the *forward* Hessian 
in the case that tableflag>0, but this has been changed since then).

.\"                *=======================================*
.\"================|    Known problems with the command    |===================
.\"                *=======================================*
.SH BUGS

Please DO NOT "fix bugs" locally without telling the author about them.

The spline interpolation may cause artificial low velocity zones if you feed the program
a model with a very steep gradient. The resulting shadow zone(s) may be very annoying,
and triplications will yield erratic behaviour of geometrical spreading.
This is probably not what you intend and can be avoided by replacing the gradient by a discontinuity,
or (but this is trickier) adding layers to force the spline to behave differently. 
The ellipticity correction is, of course, a linearized correction with all the shortcomings of
that. It does not very well handle large changes in seismic velocity in the crust, certainly not
if the crustal discontinuities are displaced by ellipticity beyond the level of the next discontinuity.
Let me know if you find a way to improve on such linearizations. Though the error may be several tenths
of a second it is a very long wavelenght error, which should not introduce artefacts in tomographic
interpretations beyond a contribution to angular order 2 terms.
Christophe Zaroli noted problems with defining ocean surface reflections like pwP in case
the background model is continental like AK135 or IASP91, even when CRUST2.0 defines an ocean
at the reflection point. This has not yet been made possible.

.\"                *=======================================*
.\"================|              Diagnostics              |===================
.\"                *=======================================*
.SH DIAGNOSTICS
Various warnings and error messages are issued to the standard error device 
by the progran.
.\"                *=======================================*
.\"================|                Examples               |===================
.\"                *=======================================*
.SH EXAMPLES
Examples of ray files:
.br

P wave
.br
P
.br
6371 1 1
.br
3482 2 1
.br
6371 5 0
.br

PcS wave (note arbitrary source depth)
.br
PcS
.br
5971 1 1
.br
3482 3 2
.br
6371 5 0
.br

pPKP wave
.br
pPKP
.br
6271 0 1
.br
6371 3 1
.br
3482 4 1
.br
1221.5 2 1
.br
3482 4 1
.br
6371 5 0
.br
.SH SOURCE CODE
The source code is located in /usr/local/soft_sismo/GLOBALSEIS/src; the
script to recompile is in /usr/local/soft_sismo/GLOBALSEIS/compile_scripts.



.\"                *=======================================*
.\"================|Other relevant commands, files, manuals|===================
.\"                *=======================================*
.SH SEE ALSO

For the theory behind the program see:
.br
Dahlen, F.A., S.-H. Hung and G. Nolet, Frechet kernels for finite-frequency 
traveltimes - I. Theory, Geophys. J. Int., 141, 157-174, 2000.
.br
and
.br
Nolet, G., A Breviary for Seismic Tomography, Cambridge Univ. Press, 2008.
( http://www.cup.cam.ac.uk/catalogue/email.asp?isbn=9780521882446 ) 

For test results of the software and a description of the code (***this is the 
paper you should refer to when using the code***):
.br
Tian, Y., G. Nolet, S.-H. Hung, R. Montelli and F.A. Dahlen, Dynamic ray 
tracing and travel time corrections for global seismic tomography, 
J. Comp. Physics, 226, 672-687, 2007.

The crustal corrections use CRUST2.0. For this please reference:
.br
Bassin, C., Laske, G. and Masters, G., The Current Limits of Resolution for 
Surface Wave Tomography in North America, EOS Trans AGU, 81, F897, 2000,
and/or Gabi Laske's website mahi.ucsd.edu/Gabi/rem.dir/crust/crust2.html.

A more extensive description of its predecessor, CRUST 5.1, can be found in: 
.br
Mooney, Laske and Masters, Crust 5.1: a global crustal model at 5x5 degrees, 
JGR, 103, 727-747, 1998.

Other helpful publications:
.br
Hung, S.-H. F.A. Dahlen and G. Nolet, Frechet kernels for finite-frequency 
travel times - II. Examples, Geophys. J. Int., 141, 175-203, 2000.
.br
Nolet, G., F.A. Dahlen and R.Montelli, Traveltimes and amplitudes of seismic 
waves: a re-assessment, in: A. Levander and G. Nolet (eds.), Array analysis 
of broadband seismograms, AGU Monograph Ser., 37-48, 2005.

.\"                *=======================================*
.\"================|             Author's name             |===================
.\"                *=======================================*
.SH AUTHOR
Guust Nolet and Yue Tian
