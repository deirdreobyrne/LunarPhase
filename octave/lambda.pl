#!/usr/bin/perl

use POSIX "fmod";
use Math::Trig;

open(S,"SUN.OUT");
open(M,"MOON.OUT");
$_=<S>;
$_=<S>;
$_=<S>;
$_=<M>;
$_=<M>;
$_=<M>;
open(O,">dlambda.csv");
$offset=0.0;
$is_waning=1;
$rec=0;
while ($rec < 524288) {
#while ($rec < 500) {
  $sun=<S>;
  $sun =~ s/^\s*//;
  ($st, $sl)=(split(/\s+/,$sun));
  $moon=<M>;
  $moon =~ s/^\s*//;
  ($mt, $ml)=(split(/\s+/,$moon));
  die "$st <> $mt" if ($st != $mt);
  if ($mt >= 2440587.5) {
    $t = $mt - 2440587.5;
    $t = $t * 86400 - 45.0 - 50.0*($t/36525);
    $delta = fmod($ml - $sl,360.0);
    $delta += 360.0 if ($delta < 0.0);
    if ($delta > 180.0) {
      $is_waning=1;
    } else {
      if ($is_waning) {
        $is_waning = 0;
        $offset += 360.0;
      }
    }
    printf O ("%.8f,%.10f\n",$t,(($delta+$offset)/180.0*3.14159265358979323486));
    $rec++;
  }
}

close(O);
close(S);
close(M);
