# Lunar Phases

## tl;dr - the algorithm

An algorithm for calculating the phase of the moon between 1900 and 2149 is presented below.

This algorithm actually calculates the fraction of the moon's disk that is illuminated by the sun as seen from the earth. This is equal to the fraction of the moon's diameter illuminated by the sun along every chord perpendicular to the axis of the phase.

```c
#include <math.h>

/**
 * @arg sec the number of seconds since 1970 Jan 1 00:00:00 UTC
 * @return the illuminated fraction of the moon's disk as seen from earth
 *
 * @author github.com/deirdreobyrne
 */
double getPhase(double sec) {
  double l = fmod(4.455933781445413 + sec/378902.2554950839, 2.0*M_PI);
  double d = fmod(4.847300486062353 + sec/406074.7530841359, 2.0*M_PI);
  double f = fmod(3.711741520540081 + sec/374194.9504350047, 2.0*M_PI);
  double m = fmod(6.245093325555575 + sec/5022681.132791225, 2.0*M_PI);

  return (1.0-cos( d
      + 1.098294220761289e-01 * sin(l)
      - 3.647737653400362e-02 * sin(m)
      + 2.181414668583516e-02 * sin(2.0*d-l)
      + 1.343381074079354e-02 * sin(2.0*d)
      + 3.665903202379761e-03 * sin(2.0*l)
      + 1.958316814332765e-03 * sin(d)
      - 1.240209781218944e-03 * sin(2.0*f)
      + 1.168396449662115e-03 * sin(2.0*d+l)
      + 1.106039239678838e-03 * sin(2.0*(f-d))
      + 1.066672928927640e-03 * sin(2.0*(d-l))
      + 9.823080994208807e-04 * sin(2.0*d-m-l)
      + 8.193661471249879e-04 * sin(2.0*d-m)
      - 6.220874758789685e-04 * sin(m-l)
      - 5.298160075572726e-04 * sin(m+l)
      + 2.736281166191152e-04 * sin(l-2.0*f)
      + 2.068112392924412e-04 * sin(3.0*l)
      + 1.890063962376408e-04 * sin(2.0*(2.0*d-l))
      - 1.004011788792417e-04 * sin(l+2.0*f)
      ))/2.0;
}
```

The argument is the number of seconds since 1970 Jan 1 00:00:00 UTC. The return result is the fraction of the moon's disk illuminated by the sun as seen from the earth. The accuracy is to within 0.002 during the entire period 1900-2149.

Note that the result will not be *exactly* 0.0 at new moon, 1.0 at full moon, or 0.5 at the quarters. This is (mostly) because the definition of those 4 "phases of the moon" doesn't, as such, depend on the illuminated fraction of the moon's disk.

Note also that the only way of determining if the moon is waxing or waning is to take the derivative of the function in `getPhase(double)`. Alternatively a second call to `getPhase(double)`could be made with an argument a few seconds into the future to see if the result is increasing or decreasing.

Finally, the output of `getPhase(double)`will (almost certainly) *never* be 1.0 or 0.0. This is because, at new and full moon, the moon is north or south enough that a tiny amount of it remains illuminated. Indeed if the output of `getPhase(double)`is close to 1.0, it could mean that the moon is actually eclipsed by the earth!

## Accuracy, or speed and small footprint?

If you cut off some of the terms in the equation in `getPhase(double)`you can increase the speed and reduce the footprint of the algorithm in a trade-off with accuracy. This could be useful, for instance, in a smartwatch which has only a small screen anyway in which to display the moon's phase, and which has limited resources which need to be managed.

This table gives the maximum error in the period 1900 - 2149 when some terms are omitted.

| **Omit from**          | **Accuracy** |
| ----------------------:|:------------:|
| ... sin(l)             | 0.086072     |
| ... sin(m)             | 0.036238     |
| ... sin(2.0*d-l)       | 0.019950     |
| ... sin(2.0*d)         | 0.010285     |
| ... sin(2.0*l)         | 0.005221     |
| ... sin(d)             | 0.003594     |
| ... sin(2.0*f)         | 0.002930     |
| ... sin(2.0*d+l)       | 0.002823     |
| ... sin(2.0*(f-d))     | 0.002482     |
| ... sin(2.0*(d-l))     | 0.002451     |
| ... sin(2.0*d-m-l)     | 0.002298     |
| ... sin(2.0*d-m)       | 0.002161     |
| ... sin(m-l)           | 0.002052     |
| ... sin(m+l)           | 0.001997     |
| ... sin(l-2.0*f)       | 0.001957     |
| ... sin(3.0*l)         | 0.001964     |
| ... sin(2.0*(2.0*d-l)) | 0.001954     |
| ... sin(l+2.0*f)       | 0.001953     |
| *nothing*              | 0.001951     |

So if an accuracy of 9% is sufficient, the equation can be reduced to `(1.0 - cos(d))/2.0`, and hence there is no need to calculate *f, l* or *m*.

## Derivation

Please note that I only have enough knowledge of mathematics and astronomy to be dangerous. So, I wish to present the method I used to arrive at the algorithm, as I'm quite sure it's not The Right Way To Do Things<sup>tm</sup>.

My first task was to collect a large dataset against which I could run a fitting function. My source for this dataset was the excellent [Solex](http://www.solexorb.it/) program. Using the settings in the diagram below, I gathered solar and lunar right ascension, declination, and distance data for 524,288 (2<sup>19</sup>) periods of 3 hours each starting at 1970 Jan 1 00:00:00 TDT.

![Solex settings](img/solex.png)

I ended up with two CSV files - one each for the sun and moon - containing Julian date, right ascension, declination, and distance data.

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

## Now the fun starts - fitting!

to be written

## Δt — the boogeyman of astronomical predictions

Solex works in TDT, whereas our clocks work in UTC, and the difference between the two (imaginatively called Δt) is not constant and not predictable (due to minor irregularities in the rotation of the earth caused by, amongst many other things, earthquakes). I adopted a very simple and reasonably accurate formula (for historical times, anyway) for Δt —

$$
\Delta t (seconds) = 45 + 50 { d \over 36525}
$$

where *d* is the number of days elapsed since 1970 Jan 1 00:00:00 UTC.

![delta-t](img/deltat.svg)

## Refining the result

to be written
