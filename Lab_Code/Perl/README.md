# Perl

These files have not been auto-formatted yet.

The files are similar to this autoformat setup using the `perltidy` tool:

```shell
perltidy -ce -l=0 -pt=2
```

Where

   * -ce means the "}" is on the same line as the "} else"
   * -l=0 means no maximum line length
   * -pt=2 means "function(something(complicated))" instead of "function( something(complicated) )"

The only issue is that the autoformatter makes strange decisions in some cases, such as where there are lists with comments, e.g.:

```perl
# This questionable code
  something
, thing
, another_thing   # here's another thing
, last_thing
```

is not "elegantly" formatted as

```perl
# This questionable code
  something,
  thing,
  another_thing,   # here's another thing
  last_thing
```

but instead is more like this:

```perl
# This questionable code
  something,
  thing,
  another_thing   # here's another thing
, last_thing
```

Where the comma is not able to be placed after "another_thing" (due to the trailing commenet)
