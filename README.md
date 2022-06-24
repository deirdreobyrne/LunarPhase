# Lunar Phases

## tl;dr - the algorithm

An algorithm for calculating the phase of the moon between 1900 and 2149 is presented below.

This algorithm actually calculates the fraction of the moon's disk that is illuminated by the sun as seen from the earth. This is equal to the fraction of the moon's diameter illuminated by the sun along every chord perpendicular to the axis of the phase.

```c
#include <math.h>

/**
 * @arg sec the number of seconds since 1970 Jan 1 00:00:00 UTC
 * @arg is_waxing if not null, then will contain true if the moon is waxing
 * @return the illuminated fraction of the moon's disk as seen from earth,
 *         accurate to 0.3% during the entire period 1970 - 2149
 *
 * @author github.com/deirdreobyrne
 */
double getPhase(double sec, int *is_waxing) {

  double d = fmod(4.847408287988257 + sec/406074.7465115577, 2.0*M_PI);
  double m = fmod(6.245333801867877 + sec/5022682.784840698, 2.0*M_PI);
  double l = fmod(4.456038755040014 + sec/378902.2499653011, 2.0*M_PI);
  double i = fmod(d
          +1.089809730923715e-01 * sin(l)
          -3.614132757006379e-02 * sin(m)
          +2.228248661252023e-02 * sin(2.0*d-l)
          +1.353592753655652e-02 * sin(2.0*d)
          +4.238560208195022e-03 * sin(2.0*l)
          +1.961408105275610e-03 * sin(d),2*M_PI);

  if (is_waxing) {
    *is_waxing = (i <= M_PI);
  }
  return (1.0 - cos(i))/2.0;
}
```

The argument is the number of seconds since 1970 Jan 1 00:00:00 UTC. The return result is the fraction of the moon's disk illuminated by the sun as seen from the earth. The accuracy is to within 0.003 during the entire period 1900-2149.

## Accuracy, or speed and small footprint?

If you cut off some of the terms in the equation in `getPhase(double)`you can increase the speed and reduce the footprint of the algorithm in a trade-off with accuracy. This could be useful, for instance, in a smartwatch which has only a small screen anyway in which to display the moon's phase, and which has limited resources which need to be managed.

This table gives the maximum error in the period 1900 - 2149 when some terms are omitted.

| **Omit starting from** | **Accuracy** |
| -----------------------:|:------------:|
| ... sin(l)         | 0.085961     |
| ... sin(m)         | 0.036612     |
| ... sin(2.0*d-l)   | 0.020229     |
| ... sin(2.0*d)     | 0.010333     |
| ... sin(2.0*l)     | 0.004811     |
| ... sin(d)         | 0.003269     |
| *nothing*          | 0.002875     |


## Derivation

Please note that I only have enough knowledge of mathematics and astronomy to be dangerous. So, I wish to present the method I used to arrive at the algorithm, as I'm quite sure it's not The Right Way To Do Things<sup>tm</sup>.

My first task was to collect a large dataset against which I could run a fitting function. My source for this dataset was the excellent [Solex](http://www.solexorb.it/) program. Using the settings in the diagram below, I gathered solar and lunar right ascension, declination, and distance data for 524,288 (2<sup>19</sup>) periods of 3 hours each starting at 1970 Jan 1 00:00:00 TDT.

![Solex settings](img/solex.png)

After some *vim* magic I ended up with two CSV files - one each for the sun and moon - containing Julian date, right ascension, declination, and distance data.

The next step was to import that data into [GNU Octave](https://octave.org/index.html) and convert it into a lunar phase dataset.

```matlab
format long;

printf("\nLoading data for phase angle calculation\n");
sun=csvread("SUN.csv");
moon=csvread("MOON.csv");

ang=[zeros(524288,2)];

printf("Starting...\n");

for i = 1 : 524288
  sd=deg2rad(sun(i,3));
  sr=deg2rad(sun(i,2));
  md=deg2rad(moon(i,3));
  mr=deg2rad(moon(i,2));
  psi=acos(sin(sd)*sin(md)+cos(sd)*cos(md)*cos(sr-mr));
  rs=sun(i,4)*149.59787066;
  rm=moon(i,4);
  ang(i,1)=sun(i,1)-2451545;
  ang(i,2)=psi + atan2(rm*sin(psi), rs-rm*cos(psi));
  if (mod(i,5000) == 0)
    printf("(%d / 524288)\n",i);
  endif;
end;

printf("\nDONE!\n");

csvwrite("i.csv",ang);
```

This actually creates a CSV file containing the number of days since 2000 Jan 1.5 (a standard epoch in astronomy), with a value related to the selenocentric elongation of the Earth. I had intended to calculate the latter using the following formulae (taken from Jean Meeus' excellent 1991 book *Astronomical Algorithms*)

$$
\psi = cos^{-1} ( sin (\delta_s) sin (\delta_m) + cos (\delta_s) cos (\delta_m) cos (\alpha_s - \alpha_m)) )
$$

$$
i = tan^{-1} \left( {{R_s sin (\psi)} \over {R_m - R_s cos (\psi)}} \right)
$$

where *α<sub>s</sub>, δ<sub>s</sub>* and *R<sub>s</sub>* are the sun's right ascension, declination and distance; *α<sub>m</sub>, δ<sub>m</sub>* and *R<sub>m</sub>* are the moon's right ascension, declination and distance; *ψ* is the geocentric elongation of the moon from the sun, and *i* is the selenocentric elongation of the earth from the sun. However, due to a futile attempt to make *i* continually increasing with time, I ended up using a modified version of the second formula above —

$$
i = \psi + tan^{-1} \left( {R_m sin (\psi)} \over {R_s - R_m cos(\psi)} \right)
$$

This gives an angle *i* which is the *complement* of the selenographic elongation of the earth from the sun. It is also, technically, a worse way of doing the calculation due to the difficulties in calculating the inverse cosine of values close to ±1 which lead to potential errors in *ψ*. However I believe my mistake made no substantive difference.

When *i* is the complement of the selenocentric elongation of the earth, then the phase is simply

$$
k = \left( {1 - cos(i)} \over 2 \right)
$$

## Δt — the boogeyman of astronomical predictions

Solex works in TDT, whereas our clocks work in UTC, and the difference between the two (imaginatively called Δt) is not constant and not predictable (due to minor irregularities in the rotation of the earth caused by, amongst many other things, earthquakes). I adopted a very simple and reasonably accurate formula (for historical times, anyway) for Δt —

$$
\Delta t (seconds) = 45 + 50 { d \over 36525}
$$

where *d* is the number of days elapsed since 1970 Jan 1 00:00:00 UTC.

![delta-t](img/deltat.svg)

## Fitting

In the end, I was using two different methods to fit the equation to the data.

The fastest was using the Octave `leasqr` function. This, as the name implies, gives a least square fit to the sum of the squares of the error in each of the phase values. A sample Octave script is in [octave/leasqrfit.m](octave/leasqrfit.m).

A slower way was to use the Octave `fminsearch` function. This allowed a fit which minimised the maximum error. *Note that at present (16 Jun 2022) the result given above isn't such a fit* in spite of the fact that such a fit usually gives better results. I'll get on to it... A sample Octave script is in [octave/minerrfit.m](octave/minerrfit.m).
