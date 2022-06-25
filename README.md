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

If you cut off some of the terms in the equation for `double i` in `getPhase(...)` you can increase the speed and reduce the footprint of the algorithm in a trade-off with accuracy. This could be useful, for instance, in a smartwatch which has only a small screen anyway in which to display the moon's phase, and which has limited resources which need to be managed.

This table gives the maximum error in the period 1900 - 2149 when some terms are omitted.

| **Omit starting from** | **Accuracy** |
| ----------------------:|:------------:|
| ... sin(l)             | 0.085961     |
| ... sin(m)             | 0.036612     |
| ... sin(2.0*d-l)       | 0.020229     |
| ... sin(2.0*d)         | 0.010333     |
| ... sin(2.0*l)         | 0.004811     |
| ... sin(d)             | 0.003269     |
| *nothing*              | 0.002875     |

In other words, if an accuracy of 9% is sufficient for your application, you can use `double i = fmod(d, 2*M_PI);` and get rid of the calculations for `double m = ...`  and `double l = ...`.

## Δt — the boogeyman of astronomical predictions

Astronomical predictions can only be made in a timescale that is uniform. In most cases, that means a timescale known as "TDT". Civil time is based on a different timescale, called "UTC". UTC tries to keep track of the Earth's rotation, so that it is generally around noon when the sun is at its highest. But the Earth's rotation is not uniform - things like earthquakes affect it. So, the difference between the TDT and UTC (imaginatively called Δt) is not constant and not predictable. I adopted a very simple and reasonably accurate formula (for historical times, anyway) for Δt —

$$
\Delta t (seconds) = 45 + 50 { d \over 36525}
$$

where *d* is the number of days elapsed since 1970 Jan 1 00:00:00 UTC.

![delta-t](img/deltat.svg)
