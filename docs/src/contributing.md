# Contributing

Contributing is welcomed with open arms. This helps the package become better and simply by working on it, you are increasing its popularity.

## Coding Style

One should see the [Julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/) as an introduction.
The overall style used it the [Blue style guide](https://github.com/invenia/BlueStyle).

## Exceptions To Blue

1. Instead of leading underscores for names of private functions, we keep them as is. For example `is_row_stochastic` should be called `_is_row_stochastic` under the Blue style guide. This function is simply not exported. Since Julia allows for massive code reuse, there is a chance that a developer might want to make use of our private functions (which are actually more useful than the public ones for this package).
2. We will implicitly `return nothing` in functions. The Blue style guide requires this to be explicit. This is to make code look cleaner.
3. Inline comments that do not form a full sentence can start with a small letter or a capital letter. They must not end with a full stop. This is to match with the Blue style guide and with Julia's source code itself.

## Additions To Blue

1. All struct types should inherit from an abstract type. This makes inheritance easier for other developers using our code.
2. Only have abstract types in type annotations in functions and struct types. This makes it easier for other developers using our code.
3. Try to make every element have the same character length in a matrix.
   See the following
   ```julia
   # yes
   T = [
       0.0 1.0;
       0.5 2.3;
   ]
    # yes
   T = [
       .12 .35;
       .55 .90;
   ]
    # no
   T = [
       0 1.0;
       0.5 2.3;
   ]
    # no
   T = [
       .12 .35;
       .55 .9;
   ]
   ```
   This makes reading matrices easier. If you are ever unsure, do whatever looks the best.
4. The following headings are to be used for docstrings:
   1. The methods of the function indented by 4 spaces.
   2. A brief explanation on the function.
   3. Definitions
   4. Parameters
   5. Keywords
   6. Returns
   7. Throws
   8. Notes
   9. Examples
   10. References
   Not all of these need to be present, but where applicable, they should follow this order. Headings must have a single empty line directly above them. They should come directly after a single hash (#) and space on the same line. They must have no empty lines directly below them.
