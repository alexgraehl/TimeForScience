# Quality-of-life tools #

### `trash.pl` ← `rm` ###

### terminal\_show\_16\_colors.sh ###

  * Shows the colors that a standard terminal can display

### terminal\_show\_256\_colors.pl ###

### tree\_of\_filestructure.sh ###

# Supplemental UNIX tools #

These are uniformly slower than their built-in UNIX equivalents, but usually have more features or require less setup. Speed is typically irrelevant for these operations in most cases anyway.

### `cut.pl` ←`cut` ###

Allows re-ordering of the columns in output

### `sheet.pl` ← `column` or `less -S` ###

Allows viewing of tab-delimited file in the terminal in a user-friendly manner. (Much better than `less -S`, because the columns are all lined up.)

# Database and Set Operations #

## JOIN / SELECT ##

### `join.pl` ← `join` ###

### `select.pl` ← Databse `select` operation ###

Select only lines in a file that meet certain criteria.

### `transpose.pl` ← Table/matrix transposition ###

Transposes a tab-delimited file so the bottom-left and top-right entries are exchanged (flips along the `x = -y` diagonal).

### `sets.pl` ← Set List operations ###

Very useful for converting back and forth between matrices and pairs. Used frequently in creating networks.

### `flatten.pl` & `expand.pl` ← Set List operations ###

  * `flatten.pl`:  For each row in a tab-delimited file like this:
```
Alpha   Beta   Gamma  Delta  Beta
```

It prints out:

```
Alpha  Beta
Alpha  Gamma
Alpha  Delta
Alpha  Beta
```


  * `expand.pl`: Turns a list of pairs into a first-item-major row.

For example:
```
   Alpha  Beta
   Alpha  Gamma
   Alpha  Delta
   Delta  Gamma
   Omega  Delta
```

Would be processed by expand.pl into:
```
Alpha	Beta	Gamma	Delta
Delta	Gamma
Omega	Delta
```


# Scientific Tools #

(None here for now)